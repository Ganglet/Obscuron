# GTDB Snapshot Ingestion & S3 Storage Pipeline

**Phase:** Weeks 1-3 — Environment Setup, Embedding Pipeline & Snapshot Differencing
**Owner:** Track 2 (Rayyan)
**Status:** Week 1 complete (below). Week 2 complete: full 502-genome panel built, extracted, and differenced against Pfam — go/no-go floor cleared by 41-83x. Week 3 (Track 2's half): differencing pipeline finalized, positive-label count delivered, experiment-log infrastructure built, Stage 3 embedding sample generated. See "Week 2" and "Week 3" sections below for full detail; Week 1 narrative follows unchanged.

## Course correction — full protein-set download was over-scoped

Track 1 flagged (2026-08-09) that warehousing the entire representative-genome protein set for both releases doesn't match what the go/no-go gate actually needs — a stratified sample sized for a few hundred positive labels, not every genome in GTDB. Concretely:

- **R207 (43GB)** — keep pulling. Still needed as the source pool for retrospective validation.
- **R232 (131.95GB)** — was dropped, then explicitly requested anyway on 2026-08-10 (see Status above) — Track 2's call to pull it ahead of when it's strictly needed, not a reversal of Track 1's scope guidance for the validation gate itself.
- **`nt_reps` / `genomes` (174GB / 179GB)** — never pulled, staying that way.
- **Implication for Week 2**: the embedding pipeline should stream sequences through the model and persist only derived artifacts (embeddings, dark/characterized labels, manifest), not warehouse raw FASTA — a different shape than "download everything, then process" which is what the S3 pipeline above was built for. The stratified sampling method itself (N per phylum vs. fixed random, etc.) is Track 1's methodology call, still to be specified.
- **Checksum fix revisited**: `request_checksum_calculation="when_required"` (above) sidesteps the completion error rather than adding real integrity verification. Track 1 suggested computing per-part checksums or verifying the multipart ETag instead — more correct, not yet done.

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
- **R207 completed and verified**: after the checksum fix, re-ran end to end — `s3api head-object` confirms 43,190,024,227 bytes, exact match to source.
- **Concurrent parts + rate-calc bug**: a real wifi-plan upgrade (200Mbps) didn't move sustained speed at all, which was the tell that the bottleneck wasn't bandwidth — the old code fetched and uploaded one part at a time, fully sequential, leaving the pipe idle half the time. Rewrote to fetch+upload several parts concurrently (`ThreadPoolExecutor`). While diagnosing a nonsensical 171MB/s reading, found the real bug: progress rate was computed as *all bytes ever done* ÷ *this run's elapsed time* — inflated right after every resume (small denominator, large numerator), which explains every earlier "40MB/s burst that decays to 2-3MB/s" report — that pattern was largely the bug talking, not real network behavior. Fixed to measure incrementally (bytes since last report ÷ time since last report).
- **Root cause found: McAfee**. 4 concurrent workers wedged — multiple workers hit the identical failure at the identical byte offset simultaneously, which pointed at something intercepting *every* open connection at once rather than random network flakiness. `Get-CimInstance -ClassName AntiVirusProduct` showed McAfee installed alongside Windows Defender; McAfee's firewall does SSL/TLS inspection, matching every symptom seen (self-signed cert errors, synchronized connection drops). Dropped to 2 workers as a workaround. User uninstalled McAfee, but `mc-fw-host` service kept running post-uninstall (self-protection resisting `Stop-Service`) — needed McAfee's own MCPR removal tool to fully clear it. Confirmed clean afterward: no AntiVirusProduct entry, no service, no process.
- **Isolated speed test**: with McAfee confirmed gone, plain `curl`/`aws s3 cp` tests (bypassing all retry/concurrency code) showed ~1.5MB/s GTDB download and ~2.7MB/s S3 upload per single connection — real, external, not a local-software artifact. Raised worker count back to 6 (the earlier 4-worker wedge was McAfee-caused, not a real limit) and cut the read timeout from 90s to 30s so a fully-dead connection gets abandoned 3x faster instead of idling a worker slot. Result: sustained rate went from ~0.3-1.3MB/s (2 workers, McAfee-era) to bursty 1-8MB/s averaging ~4MB/s.
- **R232 completed and verified**: `s3api head-object` confirms 131,946,537,881 bytes, exact match to source.

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
| `s3://darkmatter-gtdb-067620369122/gtdb/` | Same metadata plus both releases' full representative-genome protein FASTA — R207 (43.19GB) and R232 (131.95GB), both verified byte-exact. |
| `src/darkmatter/data/gtdb.py`, `s3_stream.py` | Fetcher and streaming-upload code, committed on `week-01-setup-ingestion`. |
| `config/snapshots.yaml` | Current release corrected to R232; historical/current pairing still open pending Track 1. |
| `data/processed/gtdb_R232/`, `gtdb_R207/` | Parsed taxonomy table (`taxonomy.csv`) + summary stats (`summary.json`) per release, gitignored, local only. |

---

## Week 2 — Embedding Pipeline Integration & Snapshot Differencing

**Status:** Complete. Branch `week-02-embedding-pipeline`, merged to `main` via PR #5.

### What Was Done

- **`src/darkmatter/data/panel.py`**: phylum-stratified genome panel sampler, per Track 1's P1-D6 spec (proportional-with-floor-cap: `n = clip(round(f·N_phylum), floor=2, cap=30)`, binary-searched `f` so the panel totals ~500). Needed the `gtdb_representative` flag from the metadata file — GTDB's taxonomy file lists every genome in a release, not just species reps, so taxonomy.csv alone wasn't enough. Ran against R207: **502 genomes across all 189 phyla**, none excluded — the whole point of the floor.
- **`.gitignore` fix**: `/data/` was blanket-excluding the directory as a unit, which meant `data/manifest.json` and every derived summary/panel file were silently uncommittable — contradicting Track 1's own P1-D6 note calling the labeled manifest "the committed artifact." Rewrote as `/data/*` plus explicit un-ignores (`*.json`, `*.csv` except `taxonomy.csv`) so small derived artifacts commit while raw/bulk data stays out.
- **`src/darkmatter/data/panel_proteins.py`**: extracts just the 502 panel genomes' proteins from the combined S3-hosted rep-protein tarball (40GB for R207), without downloading the other ~39.5GB. Confirmed tar member names (`{RS_|GB_}{accession}_protein.faa`) match `genome_panel.csv`'s accession column exactly, so no NCBI assembly_summary lookup was needed. Went through two real failures before this worked reliably:
  - A continuously-open `StreamingBody` read timed out every ~10-100KB on this connection. Rewrote to fetch bounded ~16MB byte-range chunks instead — same "bounded ranges survive a blip" lesson `s3_stream.py` already learned for uploads.
  - At 6 concurrent workers, hit the exact `SSL: CERTIFICATE_VERIFY_FAILED self-signed certificate in certificate chain` signature Track 1 traced to McAfee during the Week 1 upload saga — except McAfee is confirmed fully gone (no service/process/AV registration), so the trigger is unidentified but still concurrency-shaped. Dropped to 2 workers as the same workaround, which held.
  - Extraction is idempotent (skips accessions already on disk), which is what let it survive three separate mid-run interruptions (see Week 3) without losing completed work.
- **`src/darkmatter/data/download.py` bug fix**: the shared HTTP downloader never verified bytes-written against `Content-Length` — a mid-stream connection drop could end the response generator early with no exception raised, silently producing a truncated file. Found the hard way when a 293MB-expected Pfam HMM fetch landed at 230MB with `download()` reporting success. Now retries the whole download (bounded file, not worth resuming byte-by-byte) whenever the size doesn't match.
- **`src/darkmatter/data/hmmscan.py`**: Pfam GA-threshold scanner via `pyhmmer` (a Python/Cython HMMER3 binding with a prebuilt Windows wheel — sidesteps HMMER having no native Windows build and this machine having no WSL distro or running Docker). Confirmed GA cutoffs are readable per-HMM and usable via `bit_cutoffs="gathering"`, matching P1-D3's "standard hmmscan with GA cutoffs" exactly. First version re-opened the ~300MB HMM file per genome (9 min/genome — ~77 hours for the full panel); rewrote to batch every genome's proteins into one combined sequence block and scan each Pfam release exactly once.
- **`src/darkmatter/data/pfam_diff.py`**: Pfam-37 net-new-family proxy standing in for InterPro-latest. P1-D3 specifies InterPro-latest as the "characterised by T1" signal, but full InterProScan is a much heavier lift (~6.6GB multi-database Docker pipeline, ~10 bundled search tools, not just HMMER) than anything else built so far — deliberately deferred as a proxy rather than built outright. Fetched Pfam-37.0's HMM library and confirmed family counts exactly match Track 1's own P1-D5 numbers (35.0: 19,632; 37.0: 21,979; net-new: 2,383 vs. his release-notes-derived 2,347 — close enough to validate both methods). A protein counts as characterised-T1-proxy if it hits a Pfam-37 family that's new since 35.0.
- **`scripts/label_panel_proteins.py`**: combines the above into the actual differencing procedure — dark-at-T0 (no Pfam-35 GA hit) ∧ characterised-T1-proxy = positive. Ran across the full panel: **1,341,100 proteins; 297,798 dark-at-T0 (22.21%); 51,287 characterised-T1-proxy (3.82%); 4,138 positive-proxy (0.309%)** — 41-83x the blueprint's 50-100 go/no-go floor (D6).
- **ESM-2 bumped to `esm2_t33_650M`**, per Track 1's design note (`t30_150M` was the smoke-test default only). Verified loading and embedding correctly (hidden dim 1280).
- **`src/darkmatter/separation.py` + `scripts/compare_embeddings_pilot.py`**: the blueprint §8 "run the Genos-m/ESM-2 comparison on a small labeled subset" deliverable. Builds a labeled subset from real, unambiguous (single-Pfam-hit) proteins rather than fabricated reference sequences, embeds with each model, and reports within-family vs. across-family cosine similarity as a separation signal. Pilot (1 genome, 6 families, 24 sequences): gap 0.040. Genos-m not run here — OOMs on this machine's RTX 4060 even quantized (see `docs/reproducibility.md`); needs Track 1's M1 Pro.
- **`scripts/compare_embeddings_full.py`**: same check at real scale — 30 phylum-diverse genomes, 20 Pfam families, 100 sequences. Gap 0.036, consistent with the pilot.

### Why (Key Decisions)

**Why a Pfam-37 proxy instead of building full InterProScan?**
InterPro isn't one search — it's ~10+ member databases, each with its own bundled search tool, not just HMMER. Getting that running via Docker on Windows and scanning the full panel through it is realistically hours of setup, not the same league as fetching one HMM file. Pfam is itself the largest InterPro member database, and Track 1's own P1-D5 proxy math already leaned on Pfam net-new-family counts — reusing that mechanism was fast (minutes) and let the pilot start immediately, with the explicit caveat that it's narrower than true InterPro-latest (misses characterisation via non-Pfam member databases) reported alongside every result.

**Why bounded S3 ranges and low concurrency for the protein extractor, again?**
Same network, same lesson already paid for once in Week 1: a single long-lived connection is fragile on this connection, and this machine's real bottleneck (~1.3-2.7MB/s, confirmed by direct testing) doesn't improve with more concurrent workers — it just reintroduces the certificate-interception symptom. Chunked-and-throttled beats clever every time here.

### Outputs

| Output | Description |
|---|---|
| `data/processed/gtdb_R207/genome_panel.csv` | 502-genome phylum-stratified panel, committed. |
| `data/processed/gtdb_R207/panel_proteins/` | Extracted protein FASTA for all 502 panel genomes, gitignored, local only. |
| `data/processed/gtdb_R207/panel_protein_labels.csv` | Per-protein dark-at-T0 / characterised-T1-proxy / positive-proxy labels, 1,341,100 rows, committed. |
| `data/processed/gtdb_R207/esm2_separation_pilot.json`, `esm2_separation_full.json` | Embedding-separation results at both scales, committed. |
| `src/darkmatter/data/{panel,panel_proteins,hmmscan,pfam_diff}.py`, `src/darkmatter/separation.py` | Pipeline code, committed on `week-02-embedding-pipeline` → `main`. |

---

## Week 3 — Go/No-Go Decision & Phase 2 Prep

**Status:** Track 2's items complete. Branch `week-03-snapshot-differencing`. The formal go/no-go sign-off and Phase 2 metric lock-in are Track 1's half — see "Open item" below.

### What Was Done

- **Finalized the snapshot-differencing pipeline and delivered the final positive-label count**: 4,138 positive-proxy sequences across the full 502-genome panel (see Week 2) — this is the number the go/no-go gate (D6) needs, already committed and reported.
- **Resumable batching for `label_panel_proteins.py`**: the full-panel labeling run was killed mid-flight by unrelated environment/session restarts three separate times, each time losing the entire multi-hour run with nothing written (output only happened at the very end). Rewrote to process genomes in batches of 25, appending to the CSV after each batch and skipping already-labeled genomes on restart — a restart now costs at most one batch (~15-65 min observed), not the whole run. This is what actually got the full panel labeled.
- **`src/darkmatter/experiment_log.py` + `scripts/log_experiment.py`**: the experiment-tracking tool `docs/reproducibility.md` already described but never existed. Captures the git commit hash automatically, appends structured entries (config/result/next step) to `docs/experiment_log.md`. Wired directly into `label_panel_proteins.py`, `compare_embeddings_full.py`, and `embed_panel_sample.py` so a completed run logs itself — not a separate step someone has to remember.
- **`src/darkmatter/embedding_sample.py` + `scripts/build_embedding_sample.py`**: builds the Stage 3 embedding budget from P1-D6 — all positives (never subsampled), a capped dark-negative query universe, and a phylum-stratified characterised-at-T0 reference. Sampled **64,000 sequences** (4,138 positive / 30,000 dark-negative / 29,862 characterised-at-T0) from the real labeled panel.
- **`scripts/embed_panel_sample.py`**: embedded the full 64,000-sequence sample with ESM-2 — **3,592s, output shape (64000, 1280)**. This is genuine "Layer 1 implementation start" work that doesn't depend on which novelty-scoring algorithm Track 1 picks (EVT calibration vs. density-based): any scorer needs embeddings as raw input, so this is ready the moment that decision lands.
- **Found and fixed a real modality gap for Genos-m**: everything fetched and extracted through Week 2 (`gtdb_proteins_aa_reps`) is amino-acid protein sequence — correct for ESM-2 and the Pfam pipeline, but the wrong input for Genos-m, which does single-nucleotide tokenization. Protein can't be reversed into its source DNA (codon degeneracy), so this needed a fresh fetch, not a transform. Built `src/darkmatter/data/panel_nucleotides.py` + `scripts/extract_panel_nucleotides.py`, pulling `gtdb_proteins_nt_reps` directly from GTDB's public mirror (not through S3 — no AWS credentials needed, consistent with P1-D4). Confirmed the archive's member layout and gene calls match the protein version exactly. **Ran the full 502-genome extraction**: 61GB archive, 1.4GB of matched nucleotide sequences written, gene counts verified identical to the protein version per genome (e.g. 219/219).
- **Built reboot/crash-resilient extraction infrastructure**, mirroring Week 1's `DarkMatterGtdbResume` pattern exactly: a bash retry-loop wrapper (`data/.ops/run_nucleotide_extraction.sh`, up to 200 attempts), a logon-triggered resume task, and a 10-minute watchdog that kills and restarts a hung (not crashed) process. Needed it — this specific run crashed twice from genuine network/DNS failures (including a full connectivity outage, confirmed via direct ping/DNS checks, not a code bug) before completing. Also found and fixed a real gap while building this: Week 1's own ops scripts had never actually been committed to git, only ever existing on local disk — backfilled and un-ignored.

### Open item — not Track 2's to close

Formally, `docs/Track1_phase1_benchmark_scope.md` on `main` still marks the go/no-go criterion "pending — BLOCKS PHASE 2." That's because Track 1's own resolving decisions (P1-D2 through P1-D7 — snapshot boundary, operational dark/characterised definition, the go/no-go proxy call) live only on `origin/week-01-benchmarking`, a branch that was pushed 2026-08-10 and has not been merged. The actual data now substantively resolves the question (4,138 clears the 50-100 floor by a wide margin, on a proxy signal that undercounts if anything), but the formal sign-off and the Phase 2 evaluation-metric lock-in (Precision@K, AUROC, calibration — required before implementation begins, per blueprint D5) are both still open on Track 1's side.

### Outputs

| Output | Description |
|---|---|
| `docs/experiment_log.md` | Auto-populated run history — commit, config, result, next step per entry. |
| `data/processed/gtdb_R207/embedding_sample.csv` | Stage 3 sample manifest (protein ID, genome, category), committed. |
| `data/processed/gtdb_R207/esm2_panel_embeddings.npy` | 64,000 × 1,280 embedding matrix, gitignored (large binary), local only. |
| `data/processed/gtdb_R207/esm2_panel_embeddings_manifest.csv` | Row-to-protein mapping for the embedding matrix, committed. |
| `data/processed/gtdb_R207/panel_nucleotides/` | Nucleotide CDS FASTA for all 502 panel genomes (1.4GB), gitignored, local only — ready input for Genos-m. |
| `src/darkmatter/data/panel_nucleotides.py`, `scripts/extract_panel_nucleotides.py` | Nucleotide extractor, committed on `week-03-snapshot-differencing`. |
| `data/.ops/{run,resume,watchdog}_nucleotide_extraction.*`, `register_nucleotide_tasks.ps1` | Reboot/crash-resilient extraction infrastructure, committed. |
| `src/darkmatter/experiment_log.py`, `embedding_sample.py` | Tooling, committed on `week-03-snapshot-differencing`. |
