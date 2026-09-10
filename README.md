# Obscuron
**Calibrated novelty detection for microbial dark matter, validated retrospectively against real characterisation events.**

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](pyproject.toml)
[![PyTorch](https://img.shields.io/badge/PyTorch-%3E%3D2.4-ee4c2c)](pyproject.toml)
[![Status](https://img.shields.io/badge/status-all%205%20layers%20built%20%C2%B7%20manuscript%20pending-brightgreen)](docs/problems_and_decisions.md)

---

## What this is

Most sequenced microbial genes have no functional annotation, because annotation is normally assigned by similarity to something already known, and for a huge share of genes nothing similar has ever been characterised. Obscuron scores each such "dark matter" sequence with a calibrated novelty score — an extreme-value-theory (EVT) tail p-value over kNN distance in a protein-language-model embedding — and validates it *retrospectively*: freeze a reference database at an earlier snapshot, score what was dark then, and grade those scores against a later snapshot in which some of those sequences have since been characterised. The headline, full-scale result: **held-out-family AUROC of 0.962 (mean, layer-22 ESM-2 embeddings, 903 Pfam-35 families)** — the scorer reliably separates a withheld protein family from everything else, on real data, at scale.

---

> **Status: all technical work of a 4-phase capstone is complete on both tracks — the benchmark + go/no-go gate, the Layer-1 EVT scorer, the Layer-3 immune self/non-self layer, and the Phase-4 extension, which grew to cover all five blueprint layers: Layer 4 multi-signal convergence, Layer 5 coding-structure discrimination, and Layer 2 structure-aware embedding (ProstT5). The full-scale Genos-m comparison is done. Only the manuscript remains. This is research code from a final-year B.Tech capstone (MPSTME, NMIMS Hyderabad) targeting ACM-BCB 2027, not a packaged library.**

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

    subgraph EMBED["Embedding"]
        ESM2["ESM-2 650M (protein)<br/>leakage-clean, UniRef50 2021_04"]
        GENOSM["Genos-m 4.7B MoE (genome)<br/>leakage-controlled, GTDB R220"]
        PROSTT5["ProstT5 (structure-aware)<br/>ProtT5-XL / 3Di, encoder-only"]
    end

    subgraph SCORE["Scoring — the five layers"]
        L1["Layer 1: kNN + GPD tail (EVT)<br/>calibrated novelty score"]
        L2["Layer 2: ProstT5 structural<br/>embedding + kNN"]
        L3["Layer 3: V-detector negative selection<br/>calibrated self / non-self"]
        L5["Layer 5: codon-position + entropy<br/>coding-vs-noise"]
        L4["Layer 4: multi-signal convergence<br/>agreement across independent axes"]
    end

    subgraph EVAL["Evaluation"]
        AUROC["Held-out-family AUROC"]
        PATK["Precision@K + lift"]
        CAL["Calibration reliability"]
        CONV["Convergence (L1 &lt;-&gt; L3, multi-axis)"]
    end

    GTDB207 --> PANEL
    PFAM35 --> DIFF
    PANEL --> DIFF --> LABEL
    GTDB232 -.-> LABEL
    PFAM37 -.-> LABEL
    LABEL --> ESM2
    LABEL --> GENOSM
    LABEL --> PROSTT5
    LABEL --> L5
    ESM2 --> L1
    ESM2 --> L3
    GENOSM --> L1
    PROSTT5 --> L2
    ESM2 --> L4
    L2 --> AUROC
    L1 --> AUROC
    L1 --> PATK
    L1 --> CAL
    L3 --> AUROC
    L5 --> AUROC
    L1 --> CONV
    L3 --> CONV
    L4 --> CONV

    style L1 fill:#2b6cb0,color:#fff
    style L2 fill:#3182ce,color:#fff
    style L3 fill:#805ad5,color:#fff
    style L4 fill:#d69e2e,color:#000
    style L5 fill:#dd6b20,color:#fff
    style ESM2 fill:#2f855a,color:#fff
    style GENOSM fill:#c05621,color:#fff
    style PROSTT5 fill:#319795,color:#fff
    style CONV fill:#d69e2e,color:#000
```

Three embedding arms (protein, genome, structure-aware) feed five independently-designed layers; the convergence check where independent axes reconverge (Layer 4, plus the Layer-1/Layer-3 check) is what lets one signal cross-validate another without shared assumptions.

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
| Layer 2 — structure-aware embedding | ProstT5 (ProtT5-XL fine-tuned on Foldseek 3Di, encoder-only) gives a structure-informed embedding without full 3-D folding — the light ESMFold substitute — scored with the same held-out-family protocol. |
| Layer 4 — multi-signal convergence | Combines independent novelty axes (Layer-1 EVT + an embedding-free genomic-context axis + a sequence-composition axis) and flags high-confidence novelty only where they agree — the "independent lines converge" principle. |
| Layer 5 — coding-structure statistics | Codon-position base bias + k-mer entropy, tested against shuffled and Markov nulls, to confirm the dark genes carry genuine coding structure rather than being spurious/non-coding artifacts. |
| Evaluation harness | Held-out-family AUROC, retrospective Precision@K + lift, and a calibration reliability check — the three metrics frozen before any scorer ran — plus the Layer-1/Layer-3 convergence check and the Layer-4 multi-axis convergence. |

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

### 3. Genome-vs-protein comparison — done at full scale, ESM-2 clearly leads

Run on AWS (A10G) over the real 502-genome nucleotide panel, then compared *within one eval* against ESM-2 on the exact same 300 largest families / 6,907 proteins:

| Model (same 300-family eval) | Held-out-family AUROC (mean / median) |
|---|---|
| ESM-2 layer 22 (protein) | **0.988 / 0.996** |
| Genos-m layer 9 (genome) | 0.739 / 0.795 |

The genomic FM is a real signal (well above chance) but a distinctly weaker detector — and since Genos-m is the *leakage-controlled* arm (its GTDB R220 pretraining overlaps the benchmark window), it underperforms the leakage-clean protein FM even with a potential leakage advantage, so the ESM-2 headline is not a leakage artifact. A cross-model consistency also emerges: a mid-late layer beats the last layer for *both* models.

> **Honest scope:** capped to the 300 largest families for a reasonable cloud runtime, so both arms read slightly higher than the full 903-family eval (ESM-2 is 0.988 here vs the 0.962 headline); the *comparison* is matched. Details: [`docs/problems_and_decisions.md` § P2-D11](docs/problems_and_decisions.md).

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

### 5. Layer 4 multi-signal convergence — independent axes agree beyond chance

Three novelty axes of different kinds — Layer-1 EVT embedding novelty, an embedding-free genomic-context novelty (how dark a gene's on-contig neighbourhood is), and a sequence-composition novelty — combined so a gene is high-confidence novel only where they agree.

| Axis pair | Spearman | Reading |
|---|---|---|
| EVT ~ genomic-context | −0.05 | independent |
| genomic-context ~ composition | +0.01 | independent |
| EVT ~ composition | +0.26 | moderate (ESM-2 encodes some composition) |

The axes are largely independent, so convergence is genuine multi-evidence, not one signal restated. The 3-way convergent set is **108 genes = 3.4× more than chance would give** under independence, and every axis *depletes* the near-known positives (composition-top lift 0.26×) — so the convergent set is the high-confidence frontier the annotation pipeline leaves behind.

![Layer 4 multi-signal convergence](results/figures/layer4_convergence_summary.png)

Details: [`docs/problems_and_decisions.md` § P4-D1/D2/D4](docs/problems_and_decisions.md).

### 6. Layer 5 coding-structure — the dark genes are genuinely coding

Codon-position base bias + k-mer entropy on all 34,138 dark genes, tested against each gene's own shuffled and Markov-1 nulls:

| Signal | Value |
|---|---|
| Codon-position base bias (median: real / shuffle / Markov-1) | 0.072 / 0.011 / 0.011 |
| **Coding-vs-noise AUROC** (real vs shuffle / vs Markov-1) | **0.94 / 0.94** |

The dark matter carries genuine reading-frame structure, i.e. these are real ORFs and not spurious calls — a direct answer to the non-coding-artifact risk, and evidence the benchmark rests on real coding sequences.

![Layer 5 coding-structure vs nulls](results/figures/layer5_coding_structure_summary.png)

> **Honest scope:** true intergenic controls need full genome assemblies (not fetched); shuffled + Markov-1 nulls are the standard available substitutes. Details: [`docs/problems_and_decisions.md` § P4-D5](docs/problems_and_decisions.md).

### 7. Layer 2 structure-aware embedding (ProstT5) — strong, just below the sequence LM

ProstT5's structure-informed encoder, scored with the same held-out-family protocol on the 60 largest families (2,698 proteins, exact-count-matched to ESM-2):

| Model (same 60 families) | Held-out-family AUROC (mean / median) |
|---|---|
| ESM-2 layer 22 (sequence) | 0.993 / 0.998 |
| **ProstT5 centered (structure)** | **0.960 / 0.972** |
| Genos-m layer 9 (genome, 300-fam ref) | 0.739 / 0.795 |

Structure-aware embedding is a strong Pfam-family separator, just below the pure sequence LM and well above the genomic FM — sensible, since Pfam families are homology-defined so a sequence model is naturally strong; structure is complementary, not superior, for family separation.

![Layer 2 ProstT5 vs sequence and genomic arms](results/figures/layer2_prostt5_summary.png)

> **Honest scope:** scoped to 60 families because a 1.5B T5 encoder on MPS is slow (the blueprint's Layer-2 compute wall); the full 300-family run is a cloud afternoon. Details: [`docs/problems_and_decisions.md` § P4-D6](docs/problems_and_decisions.md).

---

## Honest limitations

- The T1 "characterised" signal is a Pfam-37 net-new-family proxy, not full InterProScan against InterPro-latest — narrower than the original design spec, so the true positive count is understated, not overstated.
- Genos-m ran at 300-family / 6,907-protein scale (a cloud A10G run), not the full 903-family eval, and ProstT5 (Layer 2) at 60 families — both matched against ESM-2 within-eval, but neither at the full scale of the ESM-2 headline. The full-scale runs are recoverable on a rented GPU.
- The retrospective positive set is selection-biased toward near-known genes (characterisation is homology-driven), so Precision@K measures prioritisation value, not "novelty equals characterisability" — stated explicitly, not smoothed over.
- Layer 3 (immune) underperforms Layer 1 as a standalone detector at every detector count tested up to 20,000; it's reported as a corroborating signal, not a competing one.
- The full-scale dark-query flagging result uses a stratified 14,138-of-34,138 sample of the dark-query population (all positives, 33% of dark_negatives) — a compute-time scope decision, not a methods one.
- No CI pipeline and minimal automated tests (two smoke/device tests in `tests/`) — correctness is currently established by manual runs logged in `docs/experiment_log.md`, not an automated suite.
- Layer 5's coding-vs-noise controls are shuffled and Markov-1 nulls, not true intergenic sequence — the stronger intergenic control needs full genome assemblies that weren't fetched (P1-D4).
- No wet-lab or independent biological validation of any flagged sequence — every result here is a computational prioritisation signal, not a functional claim.
- Developed and run on two personal machines (an M1 Pro laptop and an RTX 4060 laptop, 8GB VRAM), not a reproducible cloud environment — hardware-specific workarounds (fp32-only on MPS, small batch sizes) are documented but not eliminated.

---

## Repository structure

```
Obscuron/
├── src/darkmatter/          # canonical package
│   ├── data/                 # GTDB/Pfam ingestion, snapshot differencing, panel building
│   ├── embeddings/           # ESM-2, Genos-m, and ProstT5 (Layer 2) embedder backends
│   ├── scoring/               # Layer 1: kNN distance + GPD tail (EVT) novelty scorer
│   ├── immune/                # Layer 3: V-detector negative selection
│   ├── convergence/           # Layer 4: genomic-context + composition novelty, multi-axis convergence
│   ├── statistical/           # Layer 5: codon-position bias + k-mer entropy coding-structure stats
│   └── analysis/              # figure generation for the evaluation harness
├── scripts/                  # CLI entry points: fetch, embed, score, evaluate, sweep, convergence, layer2/5
├── config/                   # frozen hyperparameters (snapshots, scorer, immune, convergence, statistical, models .yaml)
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

The Genos-m comparison did use a rented cloud GPU (AWS A10G) for the 300-family run; the one remaining paid-compute item is scaling it to the full 903-family eval (and ProstT5 to full scale), each roughly a ~$2 rented-GPU afternoon.

---

## Roadmap

| Phase | Weeks | Track 1 (design / analysis) | Track 2 (implementation) |
|---|---|---|---|
| 1 — Benchmark + go/no-go gate | 1–3 | ✅ Complete | ✅ Complete |
| 2 — Layer 1 EVT scorer | 4–6 | ✅ Complete | ✅ Complete |
| 3 — Layer 3 immune layer | 7–9 | ✅ Complete | ✅ Complete |
| 4 — Extension + manuscript | 10–12 | ✅ Extension complete (Layers 2, 4, 5 built; robustness + full Genos-m comparison; lit-search re-run) · 🔄 manuscript | ✅ Figures, README, citations |

---

## Authors

- **Rayyan Mohammed** — Track 2: implementation & experimental execution. Ingestion/snapshot-differencing pipeline, embedding extraction, EVT scorer implementation, immune-layer sweep and full-scale flagging infrastructure, reproducibility tooling.
- **Angshuman Chakravertty** — Track 1: methodology, design & analysis. Snapshot boundary and provenance standards, EVT/GPD scorer design, immune self/non-self layer design and build, second-boundary robustness check, the within-eval genome-vs-protein comparison, Layer 4 multi-signal convergence, Layer 5 coding-structure discrimination, the Layer 2 ProstT5 structure-aware arm, the pre-submission literature-search re-run, and manuscript authoring.

---

## Documentation index

- [`docs/problems_and_decisions.md`](docs/problems_and_decisions.md) — full numbered decision log, every design choice and result
- [`docs/Track1_phase1_benchmark_scope.md`](docs/Track1_phase1_benchmark_scope.md) — Phase 1 working record
- [`docs/Track1_phase2_scorer_design.md`](docs/Track1_phase2_scorer_design.md) — Layer 1 EVT scorer design spec
- [`docs/Track1_phase3_immune_design.md`](docs/Track1_phase3_immune_design.md) — Layer 3 immune design spec and cross-track hand-offs
- [`docs/Track1_phase4_extension_design.md`](docs/Track1_phase4_extension_design.md) — Phase 4 extension: Layer 4 convergence (+ Layers 2 & 5) design and results
- [`docs/Track2_Phase1_Execution.md`](docs/Track2_Phase1_Execution.md) — Track 2's Phase 1 execution notes
- [`docs/Track2_Phase2_scoring_handoff.md`](docs/Track2_Phase2_scoring_handoff.md) — scorer implementation interface contract
- [`docs/experiment_log.md`](docs/experiment_log.md) — one dated entry per meaningful run
- [`docs/reproducibility.md`](docs/reproducibility.md) — hardware findings and dataset provenance

---

## License

MIT — see [`LICENSE`](LICENSE). Reference materials and external tools are credited in [`ACKNOWLEDGEMENTS.md`](ACKNOWLEDGEMENTS.md).
