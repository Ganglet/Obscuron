# Phase 4 Track 1 — Extension: Layer 4 Multi-Signal Convergence

**Phase:** 4 — Extension, Analysis & Manuscript (weeks 10–12)
**Track:** 1 — Methodology, Design & Analysis
**Status:** BUILT — 3-way (P4-D2 + P4-D4). Axes = ESM-2 EVT embedding novelty + genomic-context novelty + compositional/statistical novelty (the last local, closing blueprint Layer 5 partially; the Genos-m axis is cloud-blocked but wired via `--extra-axis`). Axes largely independent (evt~context −0.05, context~composition +0.01, evt~composition +0.26); all deplete near-known positives (context 0.68×, composition 0.26×, 3-way-convergent 108 genes 0.53×) — every independent novelty line points to the frontier. Convergence beyond chance (108 vs ~32 expected under independence).
**Branch:** `phase-3-track-1-Immune-inspired-self/non-self-discrimination` (Phase-4 work continues here until the phase-4 branch is cut)

---

## Objective

Deliver the blueprint's Phase-4 extension (§7 "a limited extension from Layer 2 or
Layer 4, if time permits") as **Layer 4 — multi-signal convergence** (§6). Combine
independent novelty axes and assign high confidence **only where they agree**, the
§4 "independent lines converge" principle that the whole framework is built on.

Layer 4 was chosen over Layer 2 (structural fusion) because Layer 2 is the one
compute-constrained path in the project (ESMFold/ProstT5 structure prediction, §9),
whereas Layer 4 reuses signals already in hand and adds a genuinely **independent**,
embedding-free axis at near-zero cost. It also directly extends the Layer-3
convergence result (P3-D6) into a full multi-axis convergence.

Blueprint references: §4 (convergence), §6 Layer 4, §7 (extension if time permits),
§10 (selection-bias limitation). Decisions: **P4-D1** (design), **P4-D2** (result),
building on **P2-D6** (EVT scorer + the Precision@K inversion) and **P1-D3** (labels).

Five deliverables, in dependency order:

1. Signal set — the independent novelty axes.
2. Genomic-context novelty — construction from contig gene-order.
3. Convergence method — agreement across axes.
4. Evaluation — independence + what convergence buys.
5. Frozen hyperparameters.

---

## 1. Signal set

**Status:** DECIDED — see P4-D1.

Two axes, deliberately of **different kinds** so their agreement is real evidence:

- **Axis A — sequence-embedding novelty (Layer 1).** The frozen ESM-2 EVT novelty
  score per dark query (`esm2_novelty_scores.csv`, `evt_novelty_score`). Distance
  from characterised sequence space.
- **Axis B — genomic-context novelty (new).** How dark a gene's genomic
  neighbourhood is — embedding-free, computed from GTDB gene-order (deliverable 2).

- [x] Two independent axes fixed (A = ESM-2 EVT, B = genomic context).
- [x] **Axis C — compositional/statistical novelty** BUILT (P4-D4): Mahalanobis of AA
  composition vs the characterised reference. Local, embedding-free, closes Layer 5
  partially. Chosen as the deliverable third axis because it is buildable now.
- [x] Genos-m per-query novelty scoped as an **optional fourth axis** — cloud-blocked
  (needs Genos-m dark-query + reference embeddings on AWS, P2-D11); **wired** via the
  generic `--extra-axis <csv>` mechanism, so a Rayyan cloud run drops straight in.

**Why not Layer 3 as axis B:** the Layer-3 negative-selection score is also
ESM-2-embedding-based (correlated +0.32 with Layer 1, P3-D6), so it is a weaker
independent axis than genomic context, which shares no input with the embedding.

---

## 2. Genomic-context novelty

**Status:** DECIDED — see P4-D1. **BUILT** (`darkmatter/convergence/multisignal.py`).

Guilt-by-association, inverted into a novelty signal. A dark gene surrounded by
**characterised** genes sits in a well-understood neighbourhood (a known operon /
pathway context) → more tractable → **low** context novelty. A dark gene embedded
among **other dark** genes sits in unexplored genomic territory → **high** context
novelty.

- **Construction:** GTDB gene-order is encoded in the protein ID
  (`<contig>_<gene-ordinal>`). Group genes by contig, order by ordinal; for each
  dark query take its ± `window` on-contig neighbours (all genes, not only scored
  queries, so the neighbourhood is the true genomic one); context novelty =
  **fraction of neighbours that are dark-at-T₀**.
- **Undefined, not zero:** genes with fewer than `min_neighbors` defined neighbours
  (contig ends, tiny contigs) are omitted — context cannot be measured, so it is not
  faked as 0.

- [x] Context = neighbourhood dark-fraction from contig gene-order (no model, no embedding).
- [x] Edge/short-contig handling = omit when under `min_neighbors`.

---

## 3. Convergence method

**Status:** DECIDED — see P4-D1.

- **Independence first:** report Spearman + top-K Jaccard between the two axes on the
  shared query set. Convergence is only meaningful multi-evidence **if the axes are
  independent** — this is the precondition, reported before anything else.
- **High-confidence convergent set:** a gene is a high-confidence novel candidate
  only if it is above the per-axis `quantile` cut on **both** axes. That is the §4
  principle made operational.

- [x] Independence (Spearman + Jaccard) is the load-bearing precondition, reported first.
- [x] Convergent set = top-quantile on BOTH axes (agreement), not a weighted sum.

---

## 4. Evaluation

**Status:** DECIDED — see P4-D1. **DONE** — P4-D2.

- **Axis independence** — is genomic context a genuinely new signal? (If Spearman ≈ 0,
  yes; convergence is real.)
- **What agreement buys** — positive-rate of the convergent set vs each single axis vs
  the base rate. Read through the **selection-bias lens (§10, P2-D6):** positives are
  the *near-known* genes, so a good novelty signal is *expected* to deplete them — the
  claim is prioritisation of the genuinely-unexplored frontier, not a positive-rate lift.
- **Honest stance:** if convergence does not raise precision, that is reported plainly
  (it cannot, against a near-known label); the contribution is the independent
  multi-evidence + the high-confidence frontier set.

- [x] Independence + convergent-set characterisation + selection-bias framing fixed.
- [x] Underperformance-on-positives is an acceptable, reportable outcome.

---

## 5. Frozen hyperparameters

**Status:** DECIDED (declared pre-result) — see P4-D1. `config/convergence.yaml`.

- **`window` = 5** neighbours each side (≤ 10 total); **`min_neighbors` = 4**.
- **`quantile` = 0.90** (top 10% on each axis) for the convergent set.
- **top-K Jaccard** at `K ∈ {50, 100, 500, 1000}` (matches P2-D5).
- **seed = 42** (frozen Phase-1 panel seed).

- [x] All values declared pre-result in `config/convergence.yaml`.

---

**Phase 4 Track 1 — Layer-4 multi-signal convergence BUILT, 3-way (P4-D1, P4-D2, P4-D4).**
`darkmatter/convergence/` (genomic-context + compositional novelty + N-axis convergence),
drivers `scripts/run_convergence.py` + `scripts/compute_composition_novelty.py`, frozen
`config/convergence.yaml` — runs on the existing Layer-1 scores + panel labels + sequences,
no GPU / no re-embed. Result: the three axes are largely independent (genomic context
orthogonal to both; evt~composition moderate +0.26) so convergence is genuine
multi-evidence, and every axis points to the unexplored frontier (all deplete near-known
positives; 108-gene triple-convergent set = 3.4× beyond chance). With Layers 1, 3, and 4
(+ Layer 5 partial via the composition axis) built + the benchmark + the Genos-m
comparison + the lit-search re-run, ALL TECHNICAL WORK IS COMPLETE. Optional/cloud-gated:
Genos-m as a 4th axis (wired via `--extra-axis`), fuller Layer 5, Layer 2 structural fusion.
