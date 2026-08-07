# Phase 1 Track 1 — Retrospective Benchmark & Scope Gate

**Phase:** 1 — Benchmark Construction and Embedding Pipeline (time-boxed, weeks 1–3)
**Track:** 1 — Methodology, Design & Analysis
**Status:** In progress
**Weeks:** 1–3
**Branch:** `phase-1-track-1`

---

## Objective

Turn the blueprint's prose into a defensible, executable specification and pass the go/no-go gate before any modelling begins. Track 1 owns the *what and why*; Track 2 implements the ingestion and snapshot-differencing pipeline against this spec. Nothing is scored yet — the output of this phase is a trusted benchmark definition plus a data-backed decision on whether the primary method is viable.

Five deliverables, in dependency order:

1. Temporal snapshot boundary (Track 1 defines → Track 2 differences).
2. Embedding-extraction protocol (Track 1 specifies → Track 2 runs).
3. Positive-label size estimate + **go/no-go gate** (blocks Phase 2).
4. Systematic literature search (novelty positioning).
5. Documentation & reproducibility standards.

Blueprint references: §2.5 (retrospective validation), §6 Layer 1, §7 (recommended approach, metrics, go/no-go), §8 (work division). Decisions: D2, D3, D4, D5, D6, D7, D10.

---

## 1. Temporal snapshot boundary

**Status:** pending

Define T₀ (historical snapshot) and T₁ (current snapshot) from GTDB release history (D7), with written justification.

Design tension to resolve with numbers, not intuition:
- **Interval too short** → too few T₀-dark sequences characterised by T₁ → positive-label set below the go/no-go floor (D6).
- **Interval too long** → embedding models / tooling mismatch the era; provenance and taxonomy churn (GTDB reclassifications) muddy what "was dark matter" means at T₀.

To decide:
- [ ] Enumerate candidate GTDB release pairs (version + date) available for both snapshots.
- [ ] For each pair, project the positive-label count (feeds §3).
- [ ] Fix T₀/T₁ and record as `P1-D#` with the interval rationale.
- [ ] Define precisely what "dark matter at T₀" means operationally (no hit above threshold X against reference set Y) — this definition is the contract Track 2 codes.

**Cross-track hand-off:** Track 2 needs the exact snapshot versions, the reference set, and the similarity threshold/tool that defines "no detectable similarity" to build the snapshot-differencing pipeline.

---

## 2. Embedding-extraction protocol

**Status:** pending

Specify how a sequence becomes a vector — identically for Genos-m and ESM-2 so the comparison (D3) is fair.

To specify:
- [ ] **Model selection & versions** — Genos-m checkpoint (verify against source preprint); ESM-2 size (e.g. `esm2_t33_650M` vs larger) chosen for the 8GB/MPS budget (D9).
- [ ] **Input unit** — nucleotide contig/gene (Genos-m) vs translated protein (ESM-2). Note the asymmetry: Genos-m can see genomic context, ESM-2 sees only the translated ORF. Record the translation/ORF-calling step for the ESM-2 arm.
- [ ] **Layer / pooling strategy** — which hidden layer; mean-pool vs [CLS]/EOS vs per-residue→pooled. Keep it identical across models where architecture allows.
- [ ] **Batching & precision** — batch size, max sequence length / chunking, fp16/fp32 on MPS vs CUDA.
- [ ] **Normalisation** — L2-normalise embeddings before distance computation? (affects kNN/EVT downstream). Decide and justify.

**Cross-track hand-off:** Track 2 implements extraction to this spec and returns embedding matrices for the labelled subset so Track 1 can compare category separation.

---

## 3. Positive-label size estimate + go/no-go gate

**Status:** pending — **BLOCKS PHASE 2**

Estimate the retrospective positive-label set: sequences that were dark matter at T₀ and have been characterised by T₁.

- [ ] Produce the projected positive count for the chosen snapshot pair (§1).
- [ ] Compare against the floor: **≥ 50–100** for a stable AUROC (D6).
- [ ] Apply the go/no-go decision (D6):
  - Below floor → widen the T₀–T₁ interval, or broaden references (GTDB + Pfam/UniProt).
  - Still infeasible → promote **Layer 5** (coding-vs-noise statistics, no positive set required) to primary; re-scope Phase 2.
- [ ] **Selection-bias note (mandatory, D2):** characterised-since sequences are enriched for the *near-known*, so the positive set skews toward moderate novelty, not maximal novelty. State what Precision@K therefore actually rewards. This interpretation is a Track 1 deliverable, not a caveat to bury.

---

## 4. Systematic literature search

**Status:** pending

Confirm the novelty positioning (D4) is real, not assumed. Search PubMed, bioRxiv, arXiv, Semantic Scholar for the intersection of {open-set recognition, novelty/anomaly detection, extreme-value calibration} × {genomic/protein foundation-model embeddings, metagenomics}.

- [ ] Run the search; record queries + dated hits.
- [ ] Confirm no prior work applies a formal open-set / EVT-calibrated novelty score to a microbial *genomic* foundation model.
- [ ] Note the closest prior art and how Obscuron differs (feeds Related Work).
- [ ] Re-run immediately before manuscript submission (Phase 4) to catch work published during the project.

---

## 5. Documentation & reproducibility standards

**Status:** in progress (this scaffold)

- [x] Repo structure, README, LICENSE, `.gitignore`, `ACKNOWLEDGEMENTS.md`, decisions log.
- [ ] Dataset manifest schema — for every sequence: id, source snapshot, dark-at-T₀ flag, characterised-by-T₁ flag + label, provenance.
- [ ] Fixed seeds and version-controlled configs (`configs/`) for every run.
- [ ] Branch convention: `phase-N-track-N`; merge to `main` at phase boundaries.
- [ ] Reproducibility contract with Track 2: data built by scripts (never committed), embeddings cached under `data/`/`embeddings/` (gitignored), manifest committed.

---

## Open decisions to log as `P1-D#`

- Snapshot pair T₀/T₁ and interval rationale (§1).
- Operational definition of "dark matter at T₀" (§1).
- Embedding layer/pooling/normalisation choices (§2).
- Go/no-go outcome and any interval/reference broadening (§3).
