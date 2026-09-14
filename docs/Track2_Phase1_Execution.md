# GTDB Snapshot Ingestion & S3 Storage Pipeline

**Phase:** 1 — Environment Setup, Embedding Pipeline & Snapshot Differencing
**Owner:** Track 2 (Rayyan)
**Branch:** `week-01-setup-ingestion` -> `week-02-embedding-pipeline` -> `week-03-snapshot-differencing`
**Status:** Complete. Full 502-genome phylum-stratified panel built, extracted, and differenced against Pfam -- go/no-go floor cleared by 41-83x (4,138 positives). Stage 3 embedding sample built and embedded with ESM-2. See P1-D8, P1-D9, P1-D10 in `docs/problems_and_decisions.md` for the decision-log entries this doc summarizes and expands on.

## Objective

Configure the compute environment on both hardware targets (RTX 4060 laptop, M1 Pro), build the data-ingestion pipeline for the historical/current GTDB snapshot pair the retrospective validation strategy depends on, and deliver the differencing pipeline (dark-at-T0 vs characterised-T1-proxy vs positive-proxy) that Phase 2's scorer needs as its ground truth -- without waiting on Track 1's final snapshot-boundary sign-off, since the fetch code itself doesn't need to know the answer in advance.

## What Was Done

- **`src/darkmatter/data/gtdb.py`**: generic GTDB release fetcher, parameterized by release tag. Fetches taxonomy + quality/rep-genome metadata + provenance files (`MD5SUM.txt`, `FILE_DESCRIPTIONS.txt`) without pulling full sequence data.
- **Metadata filename bug fix**: GTDB's metadata file extension isn't stable across releases -- R232 serves `.tsv.gz`, R207 serves `.tar.gz`. Fetcher now tries `.tsv.gz` first and falls back to `.tar.gz` on a 404.
- **Release currency fix**: `config/snapshots.yaml` had `current.gtdb_release` pinned to R226, no longer GTDB's latest (R232 was). Bumped it; historical stays R207 pending Track 1's actual boundary decision.
- **Fetched and MD5-verified** both snapshots against their own `MD5SUM.txt` -- R232 (current, ~292MB metadata) and R207 (historical, ~80MB metadata). All 8 files match.
- **`.gitignore` bug fix**: the `data/` pattern (no leading slash) was matching `src/darkmatter/data/` as well as the top-level `data/` dir, silently excluding the entire data-fetcher module from version control since the first commit. Anchored to `/data/`.
- **S3 storage pipeline**: created bucket `darkmatter-gtdb-067620369122` (eu-north-1, public access blocked, AES256 default encryption). The GTDB representative-genome protein FASTA files that feed the embedding pipeline are too large for local disk (R232: ~123GB, R207: ~40GB -- neither drive on this machine has that much free space), so built `src/darkmatter/data/s3_stream.py` to stream HTTP -> S3 multipart upload directly, without buffering the full file locally.
- **IAM cleanup**: replaced root-account access-key usage with two scoped IAM users (`darkmatter-track1`, `darkmatter-track2`), each limited to read/write on just this one bucket.
- **`src/darkmatter/data/preprocess.py`**: parses the raw taxonomy TSVs into a structured table (one column per rank) plus summary counts. Confirmed real numbers: R232 has 901,341 genomes / 199,923 species vs. R207's 317,542 / 65,703 -- database size nearly tripled across the interval, a good early signal for the go/no-go positive-label estimate.
- **S3 transfer reliability, rebuilt three times under real failures**: the first live run of `s3_stream.py` died 8.5% into the 43GB R207 file on a plain network read timeout, with no way to continue except restarting from byte 0. Rewrote it to resume from the exact byte offset, off S3's own multipart-upload state (`list_multipart_uploads` / `list_parts`) rather than a local state file, so it survives a full process restart, not just an in-process retry -- confirmed GTDB's server supports `Range` requests before relying on it. A live run kept dying at a suspiciously consistent ~100-108MB into every attempt (a connection-duration limit somewhere on the network path, not a byte-count limit): switched from one open-ended streaming request to bounded ~50MB per-part `Range` requests, so every part negotiates a fresh connection. R207 fully downloaded (100%, 43.19GB) and then failed at `complete_multipart_upload` with a missing-per-part-checksum error -- a latent bug in `upload_part()` only surfacing once a transfer finally reached 100% for the first time; the already-uploaded parts couldn't be retroactively fixed, so the 43GB had to be re-transferred from zero after the fix (`request_checksum_calculation="when_required"`).
- **S3 bucket lifecycle rule**: added `AbortIncompleteMultipartUpload` at 7 days, so a stalled or abandoned upload doesn't sit there accruing storage charges indefinitely.
- **Power settings check**: confirmed sleep-after-idle is already `Never` on both AC and battery, and this machine has no lid sensor at all -- ruled out the laptop sleeping as a cause of the read-timeout crash above; it was a genuine network blip.
- **Survive a full restart, not just a process crash**: registered a logon-triggered scheduled task (`DarkMatterGtdbResume`) that runs `data/.ops/resume_transfer.ps1` 60 seconds after login, safe to fire on every login since it only relaunches when it finds the transfer stopped and incomplete.
- **Root cause of a 4-worker wedge found: McAfee.** Multiple concurrent workers hit the identical failure at the identical byte offset simultaneously -- pointed at something intercepting every open connection at once, not random flakiness. McAfee's firewall was doing SSL/TLS inspection (self-signed cert errors, synchronized connection drops matched exactly). Dropped to 2 workers as an interim workaround; user uninstalled McAfee (needed the MCPR removal tool, since its `mc-fw-host` service resisted `Stop-Service` post-uninstall); confirmed clean afterward.
- **Isolated speed test post-McAfee**: plain `curl`/`aws s3 cp` tests showed ~1.5MB/s GTDB download and ~2.7MB/s S3 upload per single connection -- real, external, not a local-software artifact. Raised worker count back to 6 (the earlier 4-worker wedge was McAfee-caused, not a real limit) and cut the read timeout from 90s to 30s. Also found and fixed a rate-calculation bug while diagnosing a nonsensical 171MB/s reading: progress rate was computed as *all bytes ever done* / *this run's elapsed time*, inflated right after every resume -- fixed to measure incrementally.
- **R207 and R232 both completed and verified**: `s3api head-object` confirms exact byte-for-byte match to source for both (43,190,024,227 bytes and 131,946,537,881 bytes).
- **`src/darkmatter/data/panel.py`**: phylum-stratified genome panel sampler, per Track 1's P1-D6 spec (proportional-with-floor-cap, binary-searched so the panel totals ~500). Ran against R207: **502 genomes across all 189 phyla**, none excluded.
- **`.gitignore` fix, again**: `/data/` was blanket-excluding the whole directory as a unit, silently making `data/manifest.json` and every derived summary/panel file uncommittable. Rewrote as `/data/*` plus explicit un-ignores so small derived artifacts commit while raw/bulk data stays out.
- **`src/darkmatter/data/panel_proteins.py`**: extracts just the 502 panel genomes' proteins from the combined S3-hosted rep-protein tarball (40GB for R207), without downloading the other ~39.5GB. Two real failures fixed before this worked reliably: a continuously-open `StreamingBody` timing out every ~10-100KB (fixed with bounded ~16MB byte-range chunks); and a certificate-verification failure at 6 concurrent workers matching McAfee's exact signature even though McAfee was confirmed gone (dropped to 2 workers, same workaround, held). Extraction is idempotent (skips accessions already on disk).
- **`src/darkmatter/data/download.py` bug fix**: the shared HTTP downloader never verified bytes-written against `Content-Length` -- a mid-stream drop could silently produce a truncated file with `download()` reporting success. Found via a 293MB-expected Pfam HMM fetch landing at 230MB. Now retries the whole download whenever the size doesn't match.
- **`src/darkmatter/data/hmmscan.py`**: Pfam GA-threshold scanner via `pyhmmer` (sidesteps HMMER having no native Windows build and no WSL/Docker on this machine). First version re-opened the ~300MB HMM file per genome (9 min/genome, ~77 hours for the full panel); rewrote to batch every genome's proteins into one combined sequence block and scan each Pfam release exactly once.
- **`src/darkmatter/data/pfam_diff.py`**: Pfam-37 net-new-family proxy standing in for InterPro-latest (full InterProScan is a much heavier ~6.6GB multi-database Docker pipeline, deliberately deferred as a proxy). Confirmed family counts exactly match Track 1's own P1-D5 numbers (35.0: 19,632; 37.0: 21,979; net-new: 2,383 vs. his 2,347).
- **`scripts/label_panel_proteins.py`**: the actual differencing procedure -- dark-at-T0 (no Pfam-35 GA hit) AND characterised-T1-proxy = positive. **Full panel result: 1,341,100 proteins; 297,798 dark-at-T0 (22.21%); 51,287 characterised-T1-proxy (3.82%); 4,138 positive-proxy (0.309%)** -- 41-83x the blueprint's 50-100 go/no-go floor. Rewritten mid-project to process genomes in batches of 25, appending to the CSV after each batch and skipping already-labeled genomes on restart, after the full-panel run was killed mid-flight by unrelated environment/session restarts three separate times and lost the entire multi-hour run each time (output only happened at the very end in the original version).
- **ESM-2 bumped to `esm2_t33_650M`** (the smoke-test default was `t30_150M` only). Verified loading and embedding correctly (hidden dim 1280).
- **`src/darkmatter/separation.py` + `scripts/compare_embeddings_{pilot,full}.py`**: the blueprint's "run the Genos-m/ESM-2 comparison on a small labeled subset" deliverable, built from real unambiguous (single-Pfam-hit) proteins rather than fabricated reference sequences. Pilot (1 genome, 6 families, 24 sequences): within-family vs. across-family cosine gap 0.040. Full (30 phylum-diverse genomes, 20 families, 100 sequences): gap 0.036, consistent with the pilot. Genos-m not run here -- OOMs on this machine's RTX 4060 even quantized; needed the M1 Pro.
- **`src/darkmatter/experiment_log.py` + `scripts/log_experiment.py`**: the experiment-tracking tool `docs/reproducibility.md` already described but never existed. Captures the git commit hash automatically, appends structured entries to `docs/experiment_log.md`, wired directly into the labeling/comparison/embedding scripts so a completed run logs itself.
- **`src/darkmatter/embedding_sample.py` + `scripts/build_embedding_sample.py`**: builds the Stage 3 embedding budget from P1-D6 -- all positives (never subsampled), a capped dark-negative query universe, a phylum-stratified characterised-at-T0 reference. Sampled **64,000 sequences** (4,138 positive / 30,000 dark-negative / 29,862 characterised-at-T0).
- **`scripts/embed_panel_sample.py`**: embedded the full 64,000-sequence sample with ESM-2 -- **3,592s, output shape (64000, 1280)**. Genuine Layer-1 groundwork that doesn't depend on which novelty-scoring algorithm Track 1 eventually picks -- any scorer needs embeddings as raw input.
- **Found and fixed a real modality gap for Genos-m**: everything fetched through this point (`gtdb_proteins_aa_reps`) is amino-acid protein sequence -- correct for ESM-2 and the Pfam pipeline, wrong input for Genos-m, which does single-nucleotide tokenization. Protein can't be reversed into its source DNA (codon degeneracy), so this needed a fresh fetch, not a transform. Built `src/darkmatter/data/panel_nucleotides.py` + `scripts/extract_panel_nucleotides.py`, pulling `gtdb_proteins_nt_reps` directly from GTDB's public mirror (no AWS credentials needed, consistent with P1-D4). Ran the full 502-genome extraction: 61GB archive, 1.4GB of matched nucleotide sequences written, gene counts verified identical to the protein version per genome (e.g. 219/219).
- **Built reboot/crash-resilient extraction infrastructure** for the nucleotide extraction, mirroring the earlier `DarkMatterGtdbResume` pattern: a bash retry-loop wrapper (up to 200 attempts), a logon-triggered resume task, and a 10-minute watchdog that kills and restarts a hung (not crashed) process. Needed it -- this run crashed twice from genuine network/DNS failures, including a full connectivity outage, before completing.
- **M1 Pro verified by Track 1**: RTX 4060 confirmed working locally (`scripts/smoke_test.py` -- ESM-2 loads and embeds; Genos-m hits the documented VRAM ceiling). Angshuman ran the M1 Pro side on his own machine: ESM-2 verified working on MPS; Genos-m deliberately left untested there rather than risk an OOM freeze on his dev machine.
- **Flagged early: full ESMFold will need cloud GPU or cluster time.** Raised per the blueprint's own instruction to flag this in Phase 1 rather than after the dataset is scoped. This phase's Genos-m findings already showed the RTX 4060's 8GB VRAM is tight even for a 4.7B-parameter embedding model with quantization -- full ESMFold structure prediction is a heavier forward pass than that. Local hardware should be assumed sufficient only for the ProstT5 lightweight substitute, not full-scale folding; any real ESMFold need in Phase 4 would require cloud GPU rental or a cluster allocation requested well before Phase 4.

## Why (Key Decisions)

**Why metadata-only first, not the full protein FASTA?**
Track 1 hadn't signed off on the historical/current snapshot pair yet. Pulling gigabytes against a boundary that might still move would've been wasted bandwidth. Metadata is enough to estimate the positive-label count for the go/no-go gate.

**Why stream directly to S3 instead of download-then-upload?**
Neither local drive has room for a 123GB protein FASTA. Streaming means bytes flow GTDB -> S3 as they arrive, no local buffering, download/upload concurrent instead of sequential.

**Why anchor `/data/` in `.gitignore`?**
An unanchored `data/` pattern matches any directory named `data` anywhere in the tree, not just the top-level one -- it was silently excluding real source code (`src/darkmatter/data/`) from version control since the first commit.

**Why scoped IAM users instead of the original access key?**
The key handed off for this was a root account key. Created per-track IAM users limited to this one bucket instead, so a leaked key can't do more than read/write GTDB snapshot data.

**Why resume from S3's multipart state instead of a local progress file?**
A local file can go stale (crash before writing it, disk not synced). S3 already tracks exactly which parts landed, for as long as the multipart upload is open -- strictly more reliable than maintaining a second copy of that fact on disk.

**Why bounded per-part ranges instead of one streaming request?**
A single connection streaming tens of GB is at the mercy of whatever kills long-lived connections on this network. Bounded ~50MB requests mean each one only needs to survive a short window.

**Why a Pfam-37 proxy instead of building full InterProScan?**
InterPro isn't one search -- it's 10+ member databases, each with its own bundled search tool. Getting that running via Docker on Windows and scanning the full panel is realistically hours of setup. Pfam is itself the largest InterPro member database, and Track 1's own P1-D5 proxy math already leaned on Pfam net-new-family counts, so reusing that mechanism was fast and let the pilot start immediately, with the narrower-than-InterPro caveat reported alongside every result.

**Why bounded S3 ranges and low concurrency for the protein extractor too?**
Same network, same lesson already paid for once: a single long-lived connection is fragile here, and the real bottleneck (~1.3-2.7MB/s, confirmed by direct testing) doesn't improve with more concurrent workers -- it just reintroduces the certificate-interception symptom.

## Outputs

| Output | Description |
|---|---|
| `data/raw/gtdb_R232/`, `data/raw/gtdb_R207/` | Fetched, MD5-verified metadata for both snapshots (gitignored, local only). |
| `data/manifest.json` | Provenance record -- source, release, fetch timestamp, file sizes for every fetch. |
| `s3://darkmatter-gtdb-067620369122/gtdb/` | Metadata plus both releases' full representative-genome protein FASTA -- R207 (43.19GB) and R232 (131.95GB), both verified byte-exact. |
| `src/darkmatter/data/{gtdb,s3_stream,panel,panel_proteins,panel_nucleotides,hmmscan,pfam_diff,preprocess,download}.py` | Fetch, streaming-upload, panel, and differencing pipeline code, committed. |
| `config/snapshots.yaml` | Current release corrected to R232; historical/current pairing left to Track 1. |
| `data/processed/gtdb_R207/genome_panel.csv` | 502-genome phylum-stratified panel, committed. |
| `data/processed/gtdb_R207/panel_proteins/`, `panel_nucleotides/` | Extracted protein and nucleotide FASTA for all 502 panel genomes, gitignored, local only. |
| `data/processed/gtdb_R207/panel_protein_labels.csv` | Per-protein dark-at-T0 / characterised-T1-proxy / positive-proxy labels, 1,341,100 rows, committed. |
| `data/processed/gtdb_R207/esm2_separation_pilot.json`, `esm2_separation_full.json` | Embedding-separation results at both scales, committed. |
| `docs/experiment_log.md` | Auto-populated run history -- commit, config, result, next step per entry. |
| `data/processed/gtdb_R207/embedding_sample.csv` | Stage 3 sample manifest (protein ID, genome, category), committed. |
| `data/processed/gtdb_R207/esm2_panel_embeddings.npy` | 64,000 x 1,280 embedding matrix, gitignored (large binary), local only. |
| `data/processed/gtdb_R207/esm2_panel_embeddings_manifest.csv` | Row-to-protein mapping for the embedding matrix, committed. |
| `data/.ops/{run,resume,watchdog}_nucleotide_extraction.*`, `register_nucleotide_tasks.ps1` | Reboot/crash-resilient extraction infrastructure, committed. |
| `src/darkmatter/experiment_log.py`, `embedding_sample.py` | Tooling, committed. |

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

**Next:** the formal go/no-go sign-off and Phase 2 evaluation-metric lock-in (Precision@K, AUROC, calibration) were Track 1's half to close -- see `docs/Track1_phase1_benchmark_scope.md`. The 4,138-positive count above (41-83x the floor, on a proxy signal that undercounts if anything) is the number that gate needed. Phase 2 scorer implementation continues in `docs/Track2_Phase2_Execution.md`.
