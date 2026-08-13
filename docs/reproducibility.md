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

## Snapshot boundary — RESOLVED (2026-08-09, Track 1)

Confirmed by **P1-D2 / P1-D3 / P1-D5** (`problems_and_decisions.md`).
**T₀ = GTDB R207 (Apr 2022) + Pfam 35.0 (Nov 2021); T₁ = GTDB R232 (Apr
2026) + InterPro-latest.** Answering the three questions this section
originally raised:

1. **Signal for "characterised since"** = Pfam/InterPro family membership,
   *not* GTDB taxonomy (GTDB supplies the sequences and a universe-growth
   number, nothing about a gene gaining function). Dark-at-T₀ = no Pfam-35
   hit at the GA threshold; characterised-by-T₁ = now hits a family new
   since T₀.
2. **Pfam→InterPro**: yes — Pfam froze at 37.0 (Jun 2024), so the *current*
   signal is **InterPro-latest**. Pfam-35 stays the historical
   dark-definition, and it aligns with ESM-2's ~2021 training cutoff,
   minimising protein-arm embedding leakage (P1-D2).
3. **Positive-label count**: proxy = +2,347 net new Pfam families 35→37 →
   ~10⁴–10⁵ positives ≫ the 50–100 floor → **GO** (P1-D5); definitive
   count = Week-2 hmmscan.

## Data storage & subsetting standard (P1-D4)

Raw reproducible public data (GTDB FASTAs) is **never warehoused** — it is
freely, permanently hosted, so the pipeline streams needed sequences from
the public mirror and persists only *derived* artifacts (embeddings,
dark/characterised labels, and a manifest recording source URLs +
checksums). The benchmark is a **principled stratified taxonomic subset** —
required for local compute (the full R207 rep set is ~10⁸ proteins) and
documented to avoid taxonomic bias. R232 proteins (application-time) and
`nt_reps`/`genomes` are deferred; any S3 staging carries a 30–60 day
lifecycle-expiry. This keeps the persistent footprint GB-scale and the S3
bill near zero.

## Embedding-model training cutoffs & per-arm leakage (P1-D7)

Load-bearing for the retrospective claim — verified, not assumed:

- **ESM-2** trained on **UniRef50 2021_04** (~late 2021). Pfam-35 was built
  on UniProt 2021_03 and T₀ = R207 (Apr 2022), so ESM-2 predates every
  T₀→T₁ characterisation event → **leakage-clean**. ESM-2 carries the
  retrospective headline.
- **Genos-m** pretrained on **GTDB R220 (released 24 Apr 2024)** + human-
  microbiome MAGs/UHGV phages (model card). R220 is *inside* the T₀→T₁
  window and is the same database the benchmark draws from → Genos-m has
  seen the benchmark sequences and 2022→Apr-2024 characterisation events
  leaked. It is self-supervised (raw DNA, no family labels) so the leak is
  weak and conservative in direction (seen sequences → lower novelty →
  harder to flag).
- **Handling:** Genos-m positives are reported **restricted to post-R220
  (Apr 2024→T₁)** for a clean number, and over the full 2022→T₁ window to
  quantify the contamination. Genos-m remains the genome-FM comparison,
  reported leakage-controlled.

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
- **Genos-m on M1 Pro — VERIFIED WORKING on MPS (2026-08-09), overturning the
  4060 result**: with the machine quiet (~11GB / 68% free), `smoke_test.py
  --skip-esm2` loaded Genos-m on `mps` (fp16, no quantization, hidden dim
  1024) and embedded toy DNA → `(3, 1024)` in 94.4s, with memory well within
  budget throughout (73% free afterward, swap ~4.6GB, no thrash). This is
  exactly the regime the 4060 could not reach: Apple unified memory has no hard
  VRAM wall, so the ~9.4GB fp16 footprint that OOM'd the 8GB 4060 — and that
  bitsandbytes could not shrink, since it cannot quantize Mixtral's
  batched-expert `Parameter` tensors — simply fits in unified memory. The
  transformers LOAD REPORT (`lm_head.weight` UNEXPECTED) is benign: base
  `MixtralModel` without the LM head; the backend mean-pools
  `last_hidden_state`. Weights are ~18GB on disk (fp32, 4 safetensors shards),
  cast to fp16 on load.
  - *Caveats*: (i) this is 3 short toy sequences — sustained throughput on real
    GTDB nucleotide data (max 8192 single-nt tokens ≈ 8kb; longer genes need
    chunking) and batch scale are untested; (ii) MPS inference on a 4.7B MoE is
    slow, so large-scale embedding may still prefer the cloud/cluster path
    (D9). But the M1 Pro is a viable Genos-m *development* machine, not
    ESM-2-only.

**Net for Track 1 development**: both backends run on this M1 Pro — ESM-2
(fast) and Genos-m (fp16 on MPS, verified). The 16GB M1 Pro clears a bar the
8GB 4060 could not, thanks to unified memory. Large-scale Genos-m throughput
and any full ESMFold remain the cloud/cluster question (blueprint §9, D9).

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
