# Illuminating Microbial Dark Matter

Structure- and novelty-aware detection of uncharacterized microbial
sequence function. Implementation of the methodology in
`Illuminating_Microbial_Dark_Matter_Blueprint.docx` — see that document
for the full research rationale; this README covers the codebase only.

## Status

Phase 1 (Weeks 1–3): embedding pipeline and data-fetch infrastructure are
built. The snapshot-boundary decision (`config/snapshots.yaml`) and the
snapshot-differencing / positive-label extraction pipeline are still
pending Track 1 sign-off — see `docs/reproducibility.md`.

## Setup

```bash
uv sync
```

Requires an NVIDIA GPU with CUDA for practical Genos-m inference (falls
back to 8-bit quantization automatically below 16GB VRAM — see
`config/models.yaml`). CPU/MPS both work for ESM-2.

## Usage

```bash
# Verify the environment and both embedding backends load on this machine
uv run python scripts/smoke_test.py

# Fetch a reference-database release (metadata/taxonomy only, not full data)
uv run python scripts/fetch_snapshot.py --source gtdb --release R226 --metadata-only
uv run python scripts/fetch_snapshot.py --source pfam --release 37.0 --metadata-only

# Embed a FASTA file with either backend
uv run python scripts/embed.py --model genos-m --fasta sequences.fasta --out out/genos_m.npy
uv run python scripts/embed.py --model esm2 --fasta sequences.fasta --out out/esm2.npy
```

## Layout

```
config/          snapshot boundaries, model/quantization policy
src/darkmatter/
  device.py        hardware detection, VRAM-aware dtype/quantization policy
  embeddings/       Genos-m + ESM-2, behind a common Embedder interface
  data/             GTDB / Pfam release fetchers, dataset manifest writer
scripts/          CLI entry points
tests/            pytest — fast (test_device.py) and slow/weight-downloading (test_embeddings_smoke.py)
docs/             reproducibility standards, experiment log
```

## Tests

```bash
uv run pytest tests/test_device.py          # fast, no GPU/weights needed
RUN_SLOW_TESTS=1 uv run pytest tests/test_embeddings_smoke.py -v   # downloads model weights
```
