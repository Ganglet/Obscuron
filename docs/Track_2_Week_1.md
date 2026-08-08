# GTDB Snapshot Ingestion & S3 Storage Pipeline

**Phase:** Week 1 — Environment Setup & Snapshot Boundary Definition
**Owner:** Track 2 (Rayyan)
**Status:** Complete on RTX 4060; M1 Pro verification still needs someone with physical access to that machine to run it (see "Open item" below). Both GTDB snapshots fetched, MD5-checked, and parsed; S3 pipeline built and transferring.

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

## Open item — M1 Pro not yet verified

RTX 4060 is confirmed working (`scripts/smoke_test.py` — ESM-2 loads and embeds; Genos-m hits the VRAM ceiling documented in `docs/reproducibility.md`). `src/darkmatter/device.py` already has an MPS code path (fp16, no quantization — bitsandbytes doesn't support MPS), but it has **not been run on real Apple Silicon hardware** — this session is on Windows, with no M1 Pro to test against. Whoever has that machine needs to run:

```bash
uv sync
uv run python scripts/smoke_test.py --skip-genos-m   # ESM-2 first; MPS + bitsandbytes can't quantize Genos-m anyway
```

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

## Outputs

| Output | Description |
|---|---|
| `data/raw/gtdb_R232/`, `data/raw/gtdb_R207/` | Fetched, MD5-verified metadata for both snapshots (gitignored, local only). |
| `data/manifest.json` | Provenance record — source, release, fetch timestamp, file sizes for every fetch. |
| `s3://darkmatter-gtdb-067620369122/gtdb/` | Same metadata plus (in progress) the full representative-genome protein FASTA for both releases. |
| `src/darkmatter/data/gtdb.py`, `s3_stream.py` | Fetcher and streaming-upload code, committed on `week-01-setup-ingestion`. |
| `config/snapshots.yaml` | Current release corrected to R232; historical/current pairing still open pending Track 1. |
| `data/processed/gtdb_R232/`, `gtdb_R207/` | Parsed taxonomy table (`taxonomy.csv`) + summary stats (`summary.json`) per release, gitignored, local only. |
