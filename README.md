# Obscuron
**Calibrated Novelty Detection for Microbial Dark Matter via Retrospective Validation**

Final-year B.Tech capstone · MPSTME, NMIMS Hyderabad. Targets ACM-BCB 2027 (fallback: IEEE BIBM).

---

## What this is

A large fraction of sequenced microbial genes match nothing in any reference database — *microbial dark matter*. Their function is unknown precisely because function is normally assigned by comparison to a known reference, and there is none. That same absence removes the ground truth needed to evaluate any method that claims to find the interesting cases.

Obscuron scores each dark-matter sequence by a **calibrated novelty score** — how far it sits outside the region of characterized sequence space in a genomic embedding — framed as open-set recognition rather than forced classification. It resolves the missing-ground-truth problem with **retrospective validation**: freeze a reference database at an earlier snapshot, score what was dark matter *then*, and grade those scores against today's more complete database, where some of those sequences have since been characterized.

The contribution is scoped as *prioritisation of sequences that merit investigation*, with an honest calibrated confidence — not functional determination, which requires wet-lab validation.

---

## Method, in one line

Embed → model the boundary of "known" → assign a calibrated novelty score (EVT / density-based) → validate retrospectively (Precision@K, held-out-family AUROC, calibration check).

- **Core (Layer 1):** calibrated novelty scoring over genomic embeddings (Genos-m primary, ESM-2 fallback *and* genome-vs-protein comparison).
- **Supporting narrative (Layer 3):** immune-inspired self/non-self discrimination over the same embeddings.
- **Extension (Layer 2/4), time permitting:** structural fusion or multi-signal convergence.

---

## Package layout and setup

Phase 1 code — embedding backends, GTDB/Pfam ingestion, snapshot fetch/differencing — lives in the `darkmatter` package under `src/darkmatter/`, built with **`uv`** (`pyproject.toml` + committed `uv.lock`). This is the canonical package; the earlier `obscuron/` conda scaffold has been removed in its favour.

```bash
uv sync

# verify the embedding backends load on this machine
uv run python scripts/smoke_test.py                  # ESM-2 + Genos-m
uv run python scripts/smoke_test.py --skip-genos-m   # ESM-2 only (fast; the M1 Pro default)

# fetch a reference-database release (metadata-only, not full sequence data)
uv run python scripts/fetch_snapshot.py --source gtdb --release R232 --metadata-only

# embed a FASTA file
uv run python scripts/embed.py --model esm2 --fasta sequences.fasta --out out/esm2.npy
```

See `docs/reproducibility.md` for hardware findings (Genos-m memory constraints on M1 Pro / RTX 4060) and dataset provenance.

---

## Data

Reference-database snapshots, sequences, and cached embeddings are **not stored in this repository** (large, and regenerable). They live under `data/`, which is gitignored, and are built by the Track 2 ingestion/snapshot-differencing pipeline. Primary reference: **GTDB** historical + current release pair; supplementary: Pfam / UniProt. See `ACKNOWLEDGEMENTS.md` for sources and citations.

---

## Structure

```
Obscuron/
├── src/darkmatter/    # canonical package: device policy, embedding backends, data pipeline
├── scripts/           # smoke test, snapshot fetch, embedding, S3 streaming
├── config/            # snapshot boundaries (snapshots.yaml), model/quantization policy (models.yaml)
├── tests/             # pytest (device, embeddings smoke)
├── data/              # gitignored — snapshots, sequences, embeddings (built locally / on S3)
├── figures/           # gitignored — output figures from evaluation
├── Paper/             # gitignored — reference PDFs (cited in ACKNOWLEDGEMENTS)
├── Documentation/     # Track 1 methodology docs + decision log
│   ├── problems_and_decisions.md              # numbered decision log (D…, P#-D…)
│   └── 01_track1_phase1_benchmark_scope.md    # Phase 1 Track 1 working record
├── docs/              # reproducibility standards, experiment log, Track 2 writeups
├── README.md
├── LICENSE
├── ACKNOWLEDGEMENTS.md
├── pyproject.toml
└── uv.lock
```

---

## Tracks

Two coordinated tracks, synchronised at each phase boundary.

- **Track 1 — Methodology, Design & Analysis** (Angshuman): snapshot boundary + provenance standards, embedding protocol + model-comparison design, metric definition, calibration methodology, result interpretation, manuscript authoring.
- **Track 2 — Implementation & Experimental Execution:** ingestion + snapshot-differencing pipeline, embedding extraction, scoring-model implementation, diagnostics, structure prediction, figures, packaging.

---

## Status

Phase 1 — Benchmark construction and scope gate (in progress).

| Phase | Weeks | Focus | Status |
|-------|-------|-------|--------|
| 1 | 1–3 *(time-boxed)* | Retrospective benchmark + embedding pipeline; **resolve go/no-go dataset gate** | **In progress** |
| 2 | 4–6 | Layer 1 calibrated novelty scorer against fixed metrics | Not started |
| 3 | 7–9 | Layer 3 immune-inspired discrimination framework | Not started |
| 4 | 10–12 | Scoped extension (if time) + manuscript | Not started |

### Phase 1 — Track 1 deliverables

- Define the temporal boundary between historical and current reference snapshots, with justification.
- Specify the embedding-extraction protocol (model selection, layer/pooling strategy, batching).
- Estimate the retrospective positive-label set size and **resolve the go/no-go criterion before Phase 2**.
- Conduct the initial systematic literature search to confirm novelty positioning.
- Establish documentation and reproducibility standards.

Working record: `Documentation/01_track1_phase1_benchmark_scope.md`.
Go/no-go floor: ≥ 50–100 characterised positive-label sequences for a stable AUROC; if short → widen the snapshot interval or broaden references; if infeasible → promote Layer 5 (coding-vs-noise statistics, which needs no positive set) to primary.

---

## License

MIT — see `LICENSE`. Reference materials and external tools are credited in `ACKNOWLEDGEMENTS.md`.
