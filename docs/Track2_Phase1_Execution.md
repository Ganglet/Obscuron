# GTDB Snapshot Ingestion & S3 Storage Pipeline

**Phase:** Week 1 — Environment Setup & Snapshot Boundary Definition
**Owner:** Track 2 (Rayyan)
**Status:** Complete on both hardware targets — RTX 4060 verified here, M1 Pro verified by Track 1 (see `docs/reproducibility.md`). Both GTDB snapshots fetched, MD5-checked, and parsed; S3 pipeline built and transferring.

## Objective

Configure the compute environment on both hardware targets (RTX 4060 laptop, M1 Pro), and build the data-ingestion pipeline for the historical/current GTDB snapshot pair the retrospective validation strategy depends on — without waiting on Track 1's final snapshot-boundary sign-off, since the fetch code itself doesn't need to know the answer in advance.

## What Was Done

- **`src/darkmatter/data/gtdb.py`**: generic GTDB release fetcher, parameterized by release tag. Fetches taxonomy + quality/rep-genome metadata + provenance files (`MD5SUM.txt`, `FILE_DESCRIPTIONS.txt`) without pulling full sequence data.
- **Metadata filename bug fix**: GTDB's metadata file extension isn't stable across releases — R232 serves `.tsv.gz`, R207 serves `.tar.gz`. Fetcher now tries `.tsv.gz` first and falls back to `.tar.gz` on a 404.
- **Release currency fix**: `config/snapshots.yaml` had `current.gtdb_release` pinned to R226, which was no longer GTDB's latest (R232 was, as of this week). Bumped it; historical stays R207 pending Track 1's actual boundary decision.
- **Fetched and MD5-verified** both snapshots against their own `MD5SUM.txt` — R232 (current, ~292MB metadata) and R207 (historical, ~80MB metadata). All 8 files match.
- **`.gitignore` bug fix**: the `data/` pattern (no leading slash) was matching `src/darkmatter/data/` as well as the top-level `data/` dir, which meant the entire data-fetcher module had never actually been committed to git since the first commit. Anchored to `/data/`.
- **S3 storage pipeline**: created bucket `darkmatter-gtdb-067620369122` (eu-north-1, public access blocked, AES256 default encryption). The GTDB representative-genome protein FASTA files that feed the embedding pipeline are too large for local disk (R232: ~123GB, R207: ~40GB — neither drive on this machine has that much free space), so built `src/darkmatter/data/s3_stream.py` to stream HTTP → S3 multipart upload directly, without buffering the full file locally.
- **IAM cleanup**: replaced root-account access-key usage with two scoped IAM users (`darkmatter-track1`, `darkmatter-track2`), each limited to read/write on just this one bucket.
- **`src/darkmatter/data/preprocess.py`**: parses the raw taxonomy TSVs into a structured table (one column per rank — domain/phylum/class/order/family/genus/species) plus summary counts. Confirmed real numbers: R232 has 901,341 genomes / 199,923 species vs. R207's 317,542 / 65,703 — database size nearly tripled across the interval, a good early signal for the go/no-go positive-label estimate.
- **S3 transfer reliability**: the first live run of `s3_stream.py` died 8.5% into the 43GB R207 file on a plain network read timeout, with no way to continue except restarting from byte 0. Rewrote it twice: first to auto-retry the whole file on failure, then properly — to resume from the exact byte offset. Resume works off S3's own multipart-upload state (`list_multipart_uploads` / `list_parts`), not a local state file, so it survives a full process restart, not just an in-process retry; confirmed GTDB's server supports `Range` requests (`Accept-Ranges: bytes`, `206` responses) before relying on it. Verified working live — a killed-and-relaunched transfer picked up at "7 parts already uploaded" instead of starting over.
- **S3 bucket lifecycle rule**: added `AbortIncompleteMultipartUpload` at 7 days, so a stalled or abandoned upload doesn't sit there accruing storage charges indefinitely; long enough not to abort a transfer that's still legitimately in progress at this connection's speed.
- **Power settings check**: confirmed sleep-after-idle is already `Never` on both AC and battery (so the machine isn't auto-sleeping mid-transfer); this machine doesn't expose a lid-close-action setting at all (no lid sensor detected), so that wasn't a factor either. The read-timeout crash above was a genuine network blip, not the laptop sleeping.
- **Survive a full restart, not just a process crash**: registered a logon-triggered scheduled task (`DarkMatterGtdbResume`) that runs `data/.ops/resume_transfer.ps1` 60 seconds after login. The script checks whether the transfer is already running or already finished and only relaunches it when it finds it stopped and incomplete — safe to fire on every login, never duplicates the process. Needed an elevated (admin) PowerShell to register; a non-elevated `schtasks /create` was denied by policy on this machine.
- **Per-part bounded range requests**: even with resume working, a live run kept dying at a suspiciously consistent ~100-108MB into every attempt — not random flakiness, but a strong sign of a connection-duration limit somewhere on the network path (router/ISP connection tracking, most likely), not a byte-count limit. The old code opened one open-ended `Range: bytes=X-` connection for the *entire remaining file* and chunked it internally, so it was always going to hit that wall. Switched to requesting each ~50MB part as its own bounded `Range: bytes=start-end`, so every part gets a fresh connection. Speed went from ~0.8MB/s-with-constant-drops to bursts of 7MB/s+ immediately after.
- **Checksum bug that broke completion after a full transfer**: R207 fully downloaded (100%, 43.19GB) and then failed at `complete_multipart_upload` with `InvalidRequest: missing checksum for part 1`. Recent botocore defaults to requiring a CRC32 checksum per part on S3 multipart uploads, but `upload_part()` here never computed/sent one — a bug latent since the very first multipart upload was created, only surfacing once a transfer finally reached 100% for the first time. The already-uploaded parts couldn't be retroactively fixed, so that upload had to be aborted and the 43GB re-transferred from zero. Fix: `Config(request_checksum_calculation="when_required")` on the S3 client, verified with a small throwaway 2-part upload before touching the real transfer again.

## M1 Pro — verified by Track 1

Closed. RTX 4060 confirmed working here (`scripts/smoke_test.py` — ESM-2 loads and embeds; Genos-m hits the VRAM ceiling documented in `docs/reproducibility.md`). Angshuman ran the M1 Pro side on his own machine (2026-08-08): ESM-2 verified working on MPS; Genos-m deliberately left untested rather than risk an OOM freeze on his dev machine — see `docs/reproducibility.md` "M1 Pro (Apple Silicon) compute environment" for the full write-up, including why unified memory makes Genos-m's fit less clear-cut than the 4060's hard VRAM wall.

## Flag — full ESMFold will need cloud GPU or cluster time

Raising this now, per the blueprint's own instruction to flag it in Phase 1 rather than after the dataset is scoped (blueprint §9, decisions log D9). This week's Genos-m findings already show the RTX 4060's 8GB VRAM is tight even for a 4.7B-parameter embedding model with quantization — full ESMFold structure prediction is a heavier forward pass than that. Local hardware (RTX 4060, and pending confirmation, M1 Pro) should be assumed sufficient only for the ProstT5 lightweight substitute and a small, high-priority candidate subset with sequence chunking, not full-scale folding. If Phase 4 ends up needing full ESMFold at any real scale, cloud GPU rental or a university cluster allocation via Dr. Maheshwari needs to be requested — that request should go out well before Phase 4, not once the extension is already underway.

## Commands

```bash
# verify the environment and both embedding backends load
uv run python scripts/smoke_test.py

# fetch a GTDB release, metadata only (taxonomy + quality flags, no sequences)
uv run python scripts/fetch_snapshot.py --source gtdb --release R232 --metadata-only
uv run python scripts/fetch_snapshot.py --source gtdb --release R207 --metadata-only

# stream the large representative-genome protein FASTA straight to S3
uv run python scripts/stream_gtdb_to_s3.py --release R232 --bucket darkmatter-gtdb-067620369122
```

## Why (Key Decisions)

**Why metadata-only first, not the full protein FASTA?**
Track 1 hadn't signed off on the historical/current snapshot pair yet (`config/snapshots.yaml` is still marked pending). Pulling gigabytes against a boundary that might still move would've been wasted bandwidth. Metadata is enough to estimate the positive-label count for the go/no-go gate.

**Why stream directly to S3 instead of download-then-upload?**
Neither local drive has room — C: is nearly full, E: has 85GB free, and the protein FASTA alone is 123GB for the current release. Streaming means bytes flow GTDB → S3 as they arrive, no local buffering, and download/upload happen concurrently instead of sequentially.

**Why anchor `/data/` in `.gitignore`?**
An unanchored `data/` pattern matches any directory named `data` anywhere in the tree, not just the top-level one. It was silently excluding `src/darkmatter/data/` — real source code — from version control since the first commit.

**Why scoped IAM users instead of the original access key?**
The key handed off for this was a root account key (unrestricted account access). Created per-track IAM users limited to this one bucket instead, so a leaked key can't do more than read/write GTDB snapshot data.

**Why CSV output for the parsed taxonomy, not parquet?**
`pandas.to_parquet` needs `pyarrow` or `fastparquet`, neither of which is a project dependency yet. CSV needed nothing new and is plenty fast at this row count (~900K rows). Worth revisiting once Week 2's differencing pipeline is built on top of it.

**Why resume from S3's multipart state instead of a local progress file?**
A local file can go stale (crash before writing it, disk not synced) and is one more thing to keep in sync with reality. S3 already tracks exactly which parts landed, for as long as the multipart upload is open — asking it directly is strictly more reliable than maintaining a second copy of that fact on disk.

**Why bounded per-part ranges instead of one streaming request?**
A single connection streaming tens of GB is at the mercy of whatever kills long-lived connections on this network. Bounded ~50MB requests mean each one only needs to survive a short window, and a fresh connection gets negotiated for every part regardless of what happened to the last one.

**Why `request_checksum_calculation="when_required"` instead of just adding checksums?**
Properly supporting checksums (compute CRC32 per part, pass it to `upload_part`, include it in `complete_multipart_upload`) is the more "correct" fix and adds real integrity verification. Opting out was faster to ship correctly and verify in isolation after already losing a full 43GB transfer to this bug once — worth revisiting if data integrity verification becomes a real requirement later.

## Outputs

| Output | Description |
|---|---|
| `data/raw/gtdb_R232/`, `data/raw/gtdb_R207/` | Fetched, MD5-verified metadata for both snapshots (gitignored, local only). |
| `data/manifest.json` | Provenance record — source, release, fetch timestamp, file sizes for every fetch. |
| `s3://darkmatter-gtdb-067620369122/gtdb/` | Same metadata plus (in progress) the full representative-genome protein FASTA for both releases. |
| `src/darkmatter/data/gtdb.py`, `s3_stream.py` | Fetcher and streaming-upload code, committed on `week-01-setup-ingestion`. |
| `config/snapshots.yaml` | Current release corrected to R232; historical/current pairing still open pending Track 1. |
| `data/processed/gtdb_R232/`, `gtdb_R207/` | Parsed taxonomy table (`taxonomy.csv`) + summary stats (`summary.json`) per release, gitignored, local only. |
