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

## Environment provenance

- `pyproject.toml` + `uv.lock` pin exact dependency versions — commit
  `uv.lock` so any teammate/reviewer can reproduce the environment with
  `uv sync`.
- `scripts/smoke_test.py` records device, VRAM, and load time for both
  embedding backends — run and note results whenever the environment or
  models.yaml quantization policy changes.

## Experiment tracking

See `docs/experiment_log.md` — one entry per meaningful run (config used,
git commit hash, result summary). Populated starting Phase 2 once the
novelty-scoring model exists to generate results worth logging.
