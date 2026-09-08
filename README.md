# Obscuron
**Calibrated novelty detection for microbial dark matter, validated retrospectively against real characterisation events.**

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](pyproject.toml)
[![PyTorch](https://img.shields.io/badge/PyTorch-%3E%3D2.4-ee4c2c)](pyproject.toml)
[![Status](https://img.shields.io/badge/status-phase%203%20complete%20%C2%B7%20phase%204%20in%20progress-yellow)](docs/problems_and_decisions.md)

---

## What this is

Most sequenced microbial genes have no functional annotation, because annotation is normally assigned by similarity to something already known, and for a huge share of genes nothing similar has ever been characterised. Obscuron scores each such "dark matter" sequence with a calibrated novelty score — an extreme-value-theory (EVT) tail p-value over kNN distance in a protein-language-model embedding — and validates it *retrospectively*: freeze a reference database at an earlier snapshot, score what was dark then, and grade those scores against a later snapshot in which some of those sequences have since been characterised. The headline, full-scale result: **held-out-family AUROC of 0.962 (mean, layer-22 ESM-2 embeddings, 903 Pfam-35 families)** — the scorer reliably separates a withheld protein family from everything else, on real data, at scale.

---

> **Status: Phases 1–3 of a 4-phase capstone are complete on both tracks (benchmark + go/no-go gate, the Layer-1 EVT scorer, and the Layer-3 immune self/non-self layer, including full-scale sweeps as of 2026-09-07). Phase 4 — folding the Genos-m and Layer-3 results into the manuscript — is in progress. This is research code from a final-year B.Tech capstone (MPSTME, NMIMS Hyderabad) targeting ACM-BCB 2027, not a packaged library.**

---

## Architecture

```mermaid
flowchart TD
    subgraph SRC["Public data sources"]
        GTDB207["GTDB R207 (T0)<br/>Apr 2022"]
        GTDB232["GTDB R232 (T1)<br/>Apr 2026"]
        PFAM35["Pfam 35.0 (T0)<br/>Nov 2021"]
        PFAM37["Pfam 37.0 proxy (T1)<br/>Jun 2024"]
    end

    subgraph INGEST["Ingestion & labeling — Track 2"]
        PANEL["502-genome phylum-stratified panel<br/>1.34M proteins scanned"]
        DIFF["Snapshot differencing"]
        LABEL["dark-at-T0 AND characterised-by-T1<br/>= 4,138 positives (0.31%)"]
    end

    subgraph EMBED["Embedding — Track 2"]
        ESM2["ESM-2 650M (protein)<br/>leakage-clean, UniRef50 2021_04"]
        GENOSM["Genos-m 4.7B MoE (genome)<br/>leakage-controlled, GTDB R220"]
    end

    subgraph SCORE["Scoring — Track 1 design, Track 2 impl"]
        L1["Layer 1: kNN + GPD tail (EVT)<br/>calibrated novelty score"]
        L3["Layer 3: V-detector negative selection<br/>calibrated self / non-self"]
    end

    subgraph EVAL["Evaluation"]
        AUROC["Held-out-family AUROC"]
        PATK["Precision@K + lift"]
        CAL["Calibration reliability"]
        CONV["Layer 1 <-> Layer 3 convergence"]
    end

    GTDB207 --> PANEL
    PFAM35 --> DIFF
    PANEL --> DIFF --> LABEL
    GTDB232 -.-> LABEL
    PFAM37 -.-> LABEL
    LABEL --> ESM2
    LABEL --> GENOSM
    ESM2 --> L1
    ESM2 --> L3
    GENOSM --> L1
    L1 --> AUROC
    L1 --> PATK
    L1 --> CAL
    L3 --> AUROC
    L1 --> CONV
    L3 --> CONV

    style L1 fill:#2b6cb0,color:#fff
    style L3 fill:#805ad5,color:#fff
    style ESM2 fill:#2f855a,color:#fff
    style GENOSM fill:#c05621,color:#fff
    style CONV fill:#d69e2e,color:#000
```

Two independently-trained embedding arms feed two independently-designed scorers (a continuous EVT distance score and a discrete immune coverage score); the convergence check where they reconverge is what lets one method cross-validate the other without shared assumptions.

---

## The problem

Annotation pipelines assign gene function by similarity search, so a gene with no similar characterised sequence gets no function — and that same absence of a known match means there's no ground truth to check whether a proposed "novelty detector" is even right. In this project's own 502-genome benchmark panel, **297,798 of 1,341,100 scanned proteins (22.21%) had zero hit against Pfam 35.0** at the family gathering threshold. Retrospective validation is the only way to manufacture a falsifiable test for a method that claims to find the interesting cases inside that 22%.

---

## How it works

| Layer / Component | What it does |
|---|---|
| Snapshot differencing | Freezes GTDB R207 + Pfam 35.0 as T0, and grades T0-dark sequences against a Pfam-37 net-new-family signal as a T1 proxy for what's been characterised since. |
| ESM-2 embedding arm | Protein-language-model embeddings (650M params, UniRef50 2021_04 cutoff), the leakage-clean headline because its pretraining predates T0. |
| Genos-m embedding arm | Genomic foundation model (4.7B MoE, nucleotide input), run as a leakage-*controlled* comparison — its GTDB R220 pretraining overlaps the benchmark window, so positives are restricted post-R220 for a clean read. |
| Layer 1 — EVT novelty scorer | kNN cosine distance to the characterised reference, calibrated with a Generalized Pareto tail fit so the output is a checked false-flag-rate p-value, not a raw distance. |
| Layer 3 — immune self/non-self | V-detector negative selection over the same embeddings: a repertoire of detectors covers "non-self" space, calibrated so the held-out-self false-flag rate stays under a target α. |
| Evaluation harness | Held-out-family AUROC, retrospective Precision@K + lift, and a calibration reliability check — the three metrics frozen before any scorer ran — plus a Layer-1/Layer-3 convergence check. |

---

## Results

### 1. The retrospective benchmark clears its go/no-go gate

| Metric | Value |
|---|---|
| Proteins scanned | 1,341,100 |
| Dark-at-T0 | 297,798 (22.21%) |
| Positive-proxy (dark → later characterised) | 4,138 (0.309%) |
| Go/no-go floor (pre-registered) | 50–100 positives |
| Margin over floor | 41–83× |

> **Honest scope:** "characterised by T1" is a Pfam-37 net-new-family proxy, not a full InterProScan run against InterPro-latest as the original design spec called for — narrower coverage, so this count is a floor, not the true number. Details: [`docs/problems_and_decisions.md` § P1-D8](docs/problems_and_decisions.md).

### 2. Layer 1 calibrated novelty scorer — 0.962 held-out-family AUROC at full scale

| Representation | AUROC (mean) | AUROC (median) | Families |
|---|---|---|---|
| Last hidden layer (raw) | 0.786 | 0.815 | 903 |
| Layer 22 (raw) | **0.962** | **0.985** | 903 |

GPD tail fit: ξ = 0.0998, β = 0.00586, KS goodness-of-fit p = 0.984.

![Held-out-family AUROC distribution](results/figures/esm2_heldout_auroc_distribution.png)
![Calibration reliability](results/figures/esm2_heldout_calibration.png)

> **Honest scope:** Precision@K lift is *below* 1 at every K tested (0.0× at K=50/100, 0.23× at K=500, 0.50× at K=1000) — this is an expected inversion, not a bug. Positives are near-known genes with intentionally *low* novelty scores, so novelty rank anti-predicts near-term characterisation. Held-out-family AUROC is the certified validation metric; Precision@K is a reported field-level finding, not a pass/fail check. Details: [`docs/problems_and_decisions.md` § P2-D6](docs/problems_and_decisions.md).

### 3. Genome-vs-protein comparison — ESM-2 leads Genos-m ~1.2× at matched scale

| Model | Separation gap (mean-centered) | Held-out AUROC (100 seqs / 20 families) |
|---|---|---|
| ESM-2 (layer 22) | 0.790 | 0.98 |
| Genos-m (layer 9) | 0.647 | 0.82 |

> **Honest scope:** this ran at a matched 100-sequence/20-family scale, not full-panel scale. Full-scale Genos-m needs more VRAM than an 8GB laptop GPU has, and the AWS free-tier plan blocks the GPU instance type that would run it — so it was deliberately deferred, not attempted and failed. Recoverable for roughly $2 on a rented GPU if ever needed. Details: [`docs/problems_and_decisions.md` § P2-D10](docs/problems_and_decisions.md).

### 4. Layer 3 immune self/non-self — calibration holds, moderate convergence with Layer 1

| Metric (α = 0.05, 5,000 detectors, PCA-50) | Value |
|---|---|
| Held-out-self flag rate (target ≤ 0.05) | 0.023 |
| Held-out-family AUROC (60 families, original build) | 0.62 mean / 0.53 median |
| Held-out-family AUROC (swept to 20,000 detectors) | 0.86 mean |
| Full-scale convergence with Layer 1 (14,138 real dark queries) | Spearman ρ = 0.485 |
| Full-scale non-self flag rate | 11.9% of dark queries — 20.1% of eventual positives vs. 8.6% of still-dark |

![Immune layer sweep summary](results/figures/immune_sweep_summary.png)

> **Honest scope:** as a standalone detector, Layer 3 stays below Layer 1 even at its best swept setting (0.86 vs. 0.962) — it was designed to *corroborate* Layer 1 via an independently-derived signal, not to beat it. The frozen production config (5,000 detectors) sits at 0.74, well below what more detectors would buy; that headroom is flagged, not yet acted on. Details: [`docs/problems_and_decisions.md` § P3-D6, P3-D7](docs/problems_and_decisions.md).

---

## Honest limitations

- The T1 "characterised" signal is a Pfam-37 net-new-family proxy, not full InterProScan against InterPro-latest — narrower than the original design spec, so the true positive count is understated, not overstated.
- Genos-m never ran at full panel scale — only a matched 100-sequence/20-family comparison exists, so the genome-vs-protein claim isn't powered at the same scale as the ESM-2 headline.
- The retrospective positive set is selection-biased toward near-known genes (characterisation is homology-driven), so Precision@K measures prioritisation value, not "novelty equals characterisability" — stated explicitly, not smoothed over.
- Layer 3 (immune) underperforms Layer 1 as a standalone detector at every detector count tested up to 20,000; it's reported as a corroborating signal, not a competing one.
- The full-scale dark-query flagging result uses a stratified 14,138-of-34,138 sample of the dark-query population (all positives, 33% of dark_negatives) — a compute-time scope decision, not a methods one.
- No CI pipeline and minimal automated tests (two smoke/device tests in `tests/`) — correctness is currently established by manual runs logged in `docs/experiment_log.md`, not an automated suite.
- No wet-lab or independent biological validation of any flagged sequence — every result here is a computational prioritisation signal, not a functional claim.
- Developed and run on two personal machines (an M1 Pro laptop and an RTX 4060 laptop, 8GB VRAM), not a reproducible cloud environment — hardware-specific workarounds (fp32-only on MPS, small batch sizes) are documented but not eliminated.

---

## Repository structure

```
Obscuron/
├── src/darkmatter/          # canonical package
│   ├── data/                 # GTDB/Pfam ingestion, snapshot differencing, panel building
│   ├── embeddings/           # ESM-2 and Genos-m embedder backends
│   ├── scoring/               # Layer 1: kNN distance + GPD tail (EVT) novelty scorer
│   ├── immune/                # Layer 3: V-detector negative selection
│   └── analysis/              # figure generation for the evaluation harness
├── scripts/                  # ~30 CLI entry points: fetch, embed, score, evaluate, sweep
├── config/                   # frozen hyperparameters (snapshots.yaml, scorer.yaml, immune.yaml, models.yaml)
├── tests/                    # pytest: device detection + embedding smoke tests
├── docs/                     # design docs, decision log, experiment log, reproducibility notes
├── results/                  # committed derived artifacts: CSVs, JSON summaries, figures
├── paper_draft/               # manuscript draft (gitignored — work in progress)
├── data/                      # gitignored — snapshots, sequences, embeddings (built locally)
├── figures/                   # gitignored — scratch evaluation output
├── pyproject.toml             # uv-managed dependencies (Python 3.11+)
└── uv.lock
```

---

## Quick start

The fastest path to seeing the core pipeline load, no cloud account and no GPU required:

```bash
git clone https://github.com/Ganglet/Obscuron.git
cd Obscuron
uv sync
uv run python scripts/smoke_test.py --skip-genos-m   # ESM-2 only — verifies the backend loads and embeds
```

---

## Reproduce everything else

```bash
# Fetch reference-database metadata (no cloud account needed)
uv run python scripts/fetch_snapshot.py --source gtdb --release R232 --metadata-only
```

```bash
# Build the benchmark panel (Track 1's frozen sampling spec, config/snapshots.yaml)
uv run python scripts/build_genome_panel.py
uv run python scripts/extract_panel_proteins.py
uv run python scripts/label_panel_proteins.py
uv run python scripts/build_embedding_sample.py
```

```bash
# Embed both arms (ESM-2 fits comfortably in 8GB VRAM; Genos-m needs ~9GB unquantized)
uv run python scripts/embed_panel_sample.py --model esm2
uv run python scripts/embed_panel_sample.py --model genos-m   # M1 Pro / 16GB+ recommended
```

```bash
# Layer 1: fit and evaluate the calibrated novelty scorer
uv run python scripts/evaluate_scorer.py --arm esm2 --config config/scorer.yaml
uv run python scripts/heldout_family_eval.py --arm esm2 --config config/scorer.yaml
```

```bash
# Layer 3: re-embed at the validated layer, build the immune detector, run the full sweep
uv run python scripts/reembed_reference_immune.py --layers 33,22 --fp32
uv run python scripts/reembed_dark_queries.py --layers 33,22 --fp32
uv run python scripts/run_immune.py
uv run python scripts/immune_sweep.py
```

---

## Data & cost hygiene

This project doesn't run live infrastructure, so the operational risk here is cloud storage cost, not uptime. **The standing rule (P1-D4): stream public data through the pipeline and persist only derived artifacts — embeddings, labels, manifests — never warehouse the raw GTDB/Pfam downloads.**

> Born from a real bill: an earlier iteration staged the full GTDB R207 protein FASTA set on S3 and racked up roughly $12 over 3 months for ~166GB of duplicated, freely re-downloadable data. If any S3 staging is used for a one-off run, apply a 30–60 day lifecycle expiry:
> ```bash
> aws s3api put-bucket-lifecycle-configuration --bucket <bucket> \
>   --lifecycle-configuration '{"Rules":[{"Expiration":{"Days":30},"Status":"Enabled"}]}'
> ```

A full-scale Genos-m run (the one piece of genuinely paid compute this project would use) is estimated at ~$2 on a rented 24GB GPU — deferred so far, not run.

---

## Roadmap

| Phase | Weeks | Track 1 (design / analysis) | Track 2 (implementation) |
|---|---|---|---|
| 1 — Benchmark + go/no-go gate | 1–3 | ✅ Complete | ✅ Complete |
| 2 — Layer 1 EVT scorer | 4–6 | ✅ Complete | ✅ Complete |
| 3 — Layer 3 immune layer | 7–9 | ✅ Complete | ✅ Complete |
| 4 — Manuscript (+ scoped extension if time) | 10–12 | 🔄 In progress | Not started |

---

## Authors

- **Rayyan Mohammed** — Track 2: implementation & experimental execution. Ingestion/snapshot-differencing pipeline, embedding extraction, EVT scorer implementation, immune-layer sweep and full-scale flagging infrastructure, reproducibility tooling.
- **Angshuman** — Track 1: methodology, design & analysis. Snapshot boundary and provenance standards, EVT/GPD scorer design, immune self/non-self layer design and initial build, manuscript authoring.

---

## Documentation index

- [`docs/problems_and_decisions.md`](docs/problems_and_decisions.md) — full numbered decision log, every design choice and result
- [`docs/Track1_phase1_benchmark_scope.md`](docs/Track1_phase1_benchmark_scope.md) — Phase 1 working record
- [`docs/Track1_phase2_scorer_design.md`](docs/Track1_phase2_scorer_design.md) — Layer 1 EVT scorer design spec
- [`docs/Track1_phase3_immune_design.md`](docs/Track1_phase3_immune_design.md) — Layer 3 immune design spec and cross-track hand-offs
- [`docs/Track2_Phase1_Execution.md`](docs/Track2_Phase1_Execution.md) — Track 2's Phase 1 execution notes
- [`docs/Track2_Phase2_scoring_handoff.md`](docs/Track2_Phase2_scoring_handoff.md) — scorer implementation interface contract
- [`docs/experiment_log.md`](docs/experiment_log.md) — one dated entry per meaningful run
- [`docs/reproducibility.md`](docs/reproducibility.md) — hardware findings and dataset provenance

---

## License

MIT — see [`LICENSE`](LICENSE). Reference materials and external tools are credited in [`ACKNOWLEDGEMENTS.md`](ACKNOWLEDGEMENTS.md).
