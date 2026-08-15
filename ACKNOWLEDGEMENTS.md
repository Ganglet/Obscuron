# Acknowledgements

Obscuron audits pretrained models and public reference databases; it trains no new
foundation model. The resources below are used under their respective licenses and
are cited here in place of vendoring them into the repository.

## Foundation models (embeddings)

**Genos-m** — genomic foundation model for human-associated microbial genomes
BGI-Research & Zhejiang Lab, bioRxiv preprint, May 2026.
Primary embedding source. **Not yet peer-reviewed** — all reported specifications
and performance figures are treated as preliminary and benchmarked independently
rather than cited from the developers (see `Documentation/problems_and_decisions.md`, D3).

**ESM-2** — protein language model
Lin et al. 2023, *Science* 379, 1123–1130. `github.com/facebookresearch/esm`
Formal fallback embedding source *and* the protein-level baseline for the
genome-vs-protein embedding comparison.

## Structure prediction (extension, Layer 2)

**ESMFold** — end-to-end structure prediction from sequence
Lin et al. 2023 (as above). Reserved for a small, high-priority candidate subset
(memory-constrained on 8GB VRAM).

**ProstT5** — structure-aware protein language model (bilingual sequence/3Di)
Heinzinger et al. 2023. Lightweight substitute for full folding.

**Foldseek** — fast structural search
van Kempen et al. 2024, *Nature Biotechnology* 42, 243–246.

## Methods referenced (no source code used)

**OpenMax** — extreme-value-theory calibration for open-set recognition
Bendale & Boult, "Towards Open Set Deep Networks," CVPR 2016.
Motivates the EVT tail-calibration of the novelty score.

**Self / non-self discrimination** — immune-inspired anomaly detection
Forrest, Perelson, Allen & Cherukuri, IEEE S&P 1994.
Basis for the Layer 3 supporting narrative.

## Reference databases

**GTDB** — Genome Taxonomy Database
Used as the historical / current snapshot pair for retrospective validation.

**Pfam / UniProt** — protein family and sequence references
Supplementary references for broadening the characterised set at the go/no-go gate.

## Key literature

- Pavlopoulos et al., "Unraveling the functional dark matter through global metagenomics," *Nature* 622, 2023 — scale of the problem.
- Chothia, "One thousand families for the molecular biologist," *Nature* 357, 1992 — structural conservation.
- "Detecting Anomalous Proteins Using Deep Representations," *NAR Genomics and Bioinformatics* 6(1), 2024 — embedding-based protein anomaly detection.

## Closest prior art (Phase-1 literature search, P1-D11)

- Ayres, Munsamy et al., "Annotating the microbial dark matter with HiFi-NN," *iScience* 2025 (PMC12148589) — closest by target; annotates dark matter to EC numbers over ESM-2 650M embeddings with a heuristic kNN confidence. Obscuron differs: open-set **novelty** (not annotation) with formal EVT calibration and a leakage-controlled retrospective benchmark.
- "DeepVirus" / "Illuminating the Virosphere's Dark Matter using Hierarchical Deep Learning," *bioRxiv* 2025 — open-set recognition + protein FM + genome context for novel viral groups. Obscuron differs: microbial **functional** dark matter (not viral lineage), EVT calibration (not hypothesis testing).
- Ma et al., "Predicting functions of uncharacterized gene products from microbial communities" (FUGAsseM), *Nature Biotechnology* 2025 — genomic-context function prediction with temporal-holdout validation. Establishes retrospective validation is *not itself novel*; Obscuron's retrospective contribution is the leakage control (P1-D7).

Full source list: blueprint §12 (`../Microbial_Dark_Matter_Blueprint.pdf`).
