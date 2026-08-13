# Reproducibility & Dataset Provenance Standards

Owned by Track 1 (methodology), maintained jointly. See blueprint §8 —
Track 1 "establishes documentation and reproducibility standards for the
dataset"; this file is where that standard lives.

## Snapshot provenance

- Every downloaded reference-database file is recorded in `data/manifest.json`
  (written automatically by `scripts/fetch_snapshot.py` via
  `darkmatter.data.manifest.write_manifest`) with source, release tag,
  fetch timestamp, and file sizes.
- The historical/current snapshot pair actually used for retrospective
  validation is the single source of truth in `config/snapshots.yaml` —
  not hardcoded anywhere else in the pipeline.
- Any change to `config/snapshots.yaml` should be a reviewed commit with a
  message explaining why the boundary moved (e.g. "widened interval,
  positive-label count was N=32, below the 50 minimum").

## Open decision: snapshot boundary

`config/snapshots.yaml` currently defaults to GTDB R207 (Apr 2022) →
R226 (Jan 2026) and Pfam 35.0 (Nov 2021) → 37.0 (Jun 2024). This is a
placeholder chosen for a reasonable (~4 year) gap, **not yet a Track 1
decision**. Before Phase 1's go/no-go gate, confirm:

1. Whether GTDB taxonomy differencing, Pfam family differencing, or both
   are the primary signal for "characterized since the historical
   snapshot."
2. Whether the Pfam→InterPro transition (Pfam has had no standalone
   release since 37.0/Jun 2024) means InterPro release notes are now a
   better source of "family added" events than Pfam directly.
3. Whether the projected positive-label count with this pair clears the
   50–100 minimum (blueprint §7) — estimate this before building the full
   differencing pipeline on top of it.

## GTDB snapshot fetch (2026-08-08)

- `config/snapshots.yaml` `current.gtdb_release` bumped R226 → R232: R226
  was no longer GTDB's latest release (confirmed against
  `data.gtdb.ecogenomic.org/releases/`, R232 published 2026-04-04). This
  is only keeping "current" current — the historical/current pairing
  itself is still Track 1's call (open decision above).
- Bug found and fixed in `src/darkmatter/data/gtdb.py`: the metadata-file
  extension is **not stable across releases** — R232 serves
  `bac120_metadata_r232.tsv.gz`, R207 serves `bac120_metadata_r207.tar.gz`.
  `fetch_gtdb_release` now tries `.tsv.gz` and falls back to `.tar.gz` on
  a 404.
- Metadata-only fetch for both R232 and R207 (taxonomy + quality/rep-genome
  metadata + MD5SUM/FILE_DESCRIPTIONS, ~370MB total) completed and MD5-verified
  against each release's `MD5SUM.txt` — all 8 files match.
- **Scale finding, relevant beyond Phase 1**: the representative-genome
  protein FASTA (`genomic_files_reps/gtdb_proteins_aa_reps_r{N}.tar.gz`) —
  what actually feeds Genos-m/ESM-2 embedding — is **~123GB on R232**, not
  the small single-file grab it might sound like. Sibling files
  (`gtdb_genomes_reps`, `gtdb_proteins_nt_reps`) are ~179GB and ~174GB.
  None of this is fetched by `metadata_only=True` or by default; a new
  `fetch_gtdb_reps_proteins()` function fetches it explicitly, on purpose,
  once storage (local disk vs. S3) and any subsetting strategy are decided.
  Budget for this before Week 2's embedding comparison.

## Genos-m hardware fit (finding from Phase 1 validation)

`scripts/smoke_test.py` was used to actually load Genos-m on the project's
RTX 4060 Laptop (8GB VRAM, 16GB system RAM) rather than assuming the
blueprint's "compatible with available hardware" claim (§9) held. It
doesn't, for a more specific reason than VRAM size alone:

- Genos-m's MoE experts are `transformers.models.mixtral.MixtralExperts` —
  batched `nn.Parameter` tensors (`gate_up_proj` shape `(32, 8192, 1024)`,
  `down_proj` shape `(32, 1024, 4096)`), not one `nn.Linear` per expert.
  bitsandbytes' `load_in_8bit`/`load_in_4bit` only replaces `nn.Linear`
  modules, so it quantizes almost none of the ~4.7B parameters (the
  experts are ~4.6B of that total) — both 8-bit and 4-bit "quantized"
  loads still tried to allocate ~6.8GB+ and hit `CUDA out of memory` on
  this GPU.
- CPU fallback (`dtype=float32`, no GPU) was also attempted and the
  process was killed without a Python traceback — consistent with an
  OOM kill rather than a code error: the model's ~9.4GB footprint plus
  Python/PyTorch overhead is tight against 16GB system RAM, and this
  machine's C: drive was nearly full at the time (see below), leaving
  little pagefile headroom to absorb the peak.

**Net effect**: GPU-resident Genos-m inference is not currently achievable
on this laptop without one of:
1. Custom quantization code that targets the batched expert `Parameter`
   tensors directly (bitsandbytes doesn't support this out of the box) —
   nontrivial engineering, not attempted.
2. More VRAM — cloud GPU rental or the university cluster allocation via
   Dr. Maheshwari, exactly the contingency blueprint §9 already names.
3. Freeing enough space on C: to give CPU inference a real chance (untested
   whether that alone would be sufficient — 9.4GB resident + overhead is
   still close to the 16GB physical RAM ceiling even with a healthy
   pagefile).

Until one of these is pursued, **ESM-2 is the practical default embedding
backend for day-to-day development on this machine** — which is exactly
the role blueprint §7 already assigns it ("If Genos-m ... becomes
unavailable at any phase boundary, ESM-2 becomes the primary embedding
source with no change to the surrounding methodology"). The ESM-2 backend
is fully verified working (`scripts/smoke_test.py --skip-genos-m`).

## M1 Pro (Apple Silicon) compute environment — verified 2026-08-08

Closes the "M1 Pro not yet verified" open item from `Track_2_Week_1.md`. Run on
the Track 1 machine (MacBook Pro, Apple M1 Pro, 10 cores, **16 GB unified
memory**, macOS) — the hardware Track 2 could not reach from the 4060.

- **Environment builds cleanly**: `uv 0.12.3` → `uv sync` on macOS with no
  manual intervention. `bitsandbytes` is correctly skipped (pyproject marker
  `platform_system != 'Darwin'`); torch 2.13.0 installs with MPS support; the
  cu121 torch index applies only to win/linux.
- **ESM-2 verified working on MPS**: `scripts/smoke_test.py --skip-genos-m`
  detects `mps`, loads `facebook/esm2_t30_150M_UR50D` (hidden dim 640, fp16,
  no quantization), embeds toy proteins → `(3, 640)`. First run 305.6s
  (dominated by the ~600MB HF download; weight-load itself <1s, cached
  thereafter). The transformers LOAD REPORT (`lm_head` UNEXPECTED, `pooler`
  MISSING) is benign — expected when an ESM-2 masked-LM checkpoint is loaded
  into `AutoModel`; the backend mean-pools `last_hidden_state` and never
  touches the pooler.
- **Genos-m on M1 Pro — not yet run definitively; a different regime from the
  4060**: on MPS `device.py` loads fp16, no quantization (~9.4GB weights).
  Unlike the 4060's hard 8GB VRAM wall, Apple unified memory has no separate
  GPU-memory ceiling — MPS can address a large fraction of the 16GB — so a fit
  is *not* foreclosed the way it is on the 4060. But 9.4GB weights +
  activations + macOS against 16GB physical is tight and may swap-thrash or be
  OOM-killed (the same ceiling that killed Track 2's CPU attempt). Left as a
  deliberate, supervised test rather than run inline, to avoid a ~9GB download
  and a possible freeze on the active development machine.

**Net for Track 1 development**: ESM-2 on MPS is the confirmed working backend
on this machine — consistent with the 4060 and with blueprint §7/D3. Genos-m
remains the cloud/cluster question (blueprint §9, D9).

## Environment provenance

- `pyproject.toml` + `uv.lock` pin exact dependency versions — commit
  `uv.lock` so any teammate/reviewer can reproduce the environment with
  `uv sync`.
- `scripts/smoke_test.py` records device, VRAM, and load time for both
  embedding backends — run and note results whenever the environment or
  models.yaml quantization policy changes.

## Experiment tracking

See `docs/experiment_log.md` — one entry per meaningful run (config used,
git commit hash, result summary, next step). `src/darkmatter/experiment_log.py`
writes entries automatically; `label_panel_proteins.py` and
`compare_embeddings_full.py` call it at the end of every run, so a
completed run is never left undocumented waiting on someone to remember
a separate logging step. Started in Phase 1 once the differencing
pipeline had a real result worth recording, not held back for Phase 2 —
`scripts/log_experiment.py` is still there for logging anything ad hoc.
