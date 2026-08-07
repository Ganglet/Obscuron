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

## Setup

```bash
conda env create -f environment.yml
conda activate obscuron
# or:  pip install -r requirements.txt
```

---

## Data

Reference-database snapshots, sequences, and cached embeddings are **not stored in this repository** (large, and regenerable). They live under `data/`, which is gitignored, and are built by the Track 2 ingestion/snapshot-differencing pipeline. Primary reference: **GTDB** historical + current release pair; supplementary: Pfam / UniProt. See `ACKNOWLEDGEMENTS.md` for sources and citations.

---

## Structure

```
Obscuron/
├── obscuron/          # source package (built phase by phase)
├── scripts/           # data build, embedding extraction, evaluation drivers
├── configs/           # version-controlled run configs (fixed seeds)
├── data/              # gitignored — snapshots, sequences, embeddings (built locally)
├── figures/           # gitignored — output figures from evaluation
├── Paper/             # gitignored — reference PDFs (cited in ACKNOWLEDGEMENTS)
├── Documentation/
│   ├── problems_and_decisions.md              # numbered decision log (D…, P#-D…)
│   └── 01_track1_phase1_benchmark_scope.md    # Phase 1 Track 1 working record
├── README.md
├── LICENSE
├── ACKNOWLEDGEMENTS.md
├── environment.yml
└── requirements.txt
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
