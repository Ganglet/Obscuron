# Phase 3 Track 1 — Immune Self/Non-Self Layer (Calibrated Negative Selection)

**Phase:** 3 — Immune-Inspired Discrimination Framework (Layer 3), weeks 7–9
**Track:** 1 — Methodology, Design & Analysis
**Status:** In progress
**Branch:** `phase-3-track-1`

---

## Objective

Integrate Layer 3 — the immune self/non-self discrimination framework (blueprint §5,
§6) — as **calibrated negative selection**: model characterised sequence space as
*self*, build a negative-selection detector repertoire that flags *non-self* (dark)
sequences with a **self-tolerance-calibrated error rate**. This is the biologically
authentic instantiation of the blueprint's Layer 3 — deliberately stronger than its
minimal "one-class kNN" suggestion — and it plays two roles: (1) the **supporting
theoretical narrative** grounding the whole method in innate-immune / negative-selection
precedent, and (2) an **independent convergence check** on the Layer-1 EVT scorer. It
**reuses the Phase-2 layer-22 embeddings at minimal extra cost** (no GPU, no re-embed).

Track 1 formalizes the framework + its evaluation + the narrative and finalizes the
config for the manuscript; Track 2 runs the parameter sweeps and scales it (blueprint §8).

Five deliverables, in dependency order:

1. Method — calibrated negative selection (the self/non-self discriminator).
2. Detector construction — V-detector negative selection in high-dimensional embedding space.
3. Self-tolerance calibration — bound the self false-flag rate at a target α.
4. Evaluation + Layer-1 convergence — held-out-family AUROC, calibration, agreement with the EVT scorer.
5. Frozen hyperparameters + immune narrative.

Blueprint references: §4 (negative selection), §5 (innate-immune precedent), §6 Layer 3,
§8 (work division). Decisions: **P3-D1…P3-D5**, building on **P2-D9** (layer-22 embeddings)
and **P2-D3/D4** (held-out-family evaluation).

---

## 1. Method — calibrated negative selection

**Status:** DECIDED — see P3-D1.

**Self = characterised-at-T₀ reference; non-self = dark-matter deviations.** A
**negative-selection detector repertoire** covers the non-self region: generate candidate
detectors in embedding space, **delete any that match self**, keep the survivors; a query
is flagged *non-self* if a surviving detector covers it. The decision threshold is
**calibrated to a self-tolerance bound** (deliverable 3), giving a binary self/non-self
call with a checked error rate plus a continuous non-self score for ranking.

This is **mechanistically distinct from Layer 1** (EVT tail-calibration on kNN distance):
Layer 1 ranks by distance to self; Layer 3 covers the *complement* with an immune
repertoire and emits a *calibrated binary* decision. So Layer 3 is additive, not a
re-narration — a different method and a different output on the same benchmark. Reuses the
**layer-22, centered ESM-2 reference embeddings** from P2-D9 (`esm2_reembed_L33_22.npz`),
so the held-out-family evaluation runs with no new embedding.

- [x] Self/non-self framing + calibrated-negative-selection method fixed (P3-D1).
- [x] Reuses P2-D9 layer-22 reference embeddings (no re-embed).
- [ ] Track 1 builds the detector + calibration + convergence analysis.

**Cross-track hand-off:** Track 1 builds `darkmatter/immune/` (detector generation,
V-detector radius, self-tolerance calibration, scoring); Track 2 runs the parameter
sweeps (§5) and scales flagging to the dark-query set.

---

## 2. Detector construction — V-detector negative selection in high-D

**Status:** DECIDED — see P3-D2.

Classic (fixed-radius) negative selection fails in 1280-D: random detectors land trivially
far from self and cover nothing useful (curse of dimensionality). Use the **V-detector
(variable-radius)** scheme instead:

- Sample a candidate centre in embedding space; set its radius to the distance to its
  **nearest self point minus a margin** (so it does not cover self); if that radius is
  below a floor, the candidate is self-covering → discard.
- Repeat until a coverage/count target is met. Survivors = the non-self repertoire.
- **High-D mitigation:** V-detector's variable radius is the primary fix; additionally
  report a variant with **PCA reduction to ~50–100 dims** (retained-variance stated) and
  compare, so the high-D behaviour is measured, not assumed.
- Geometry = cosine on the L2-normalised layer-22 embeddings (P1-D13 / P2-D9), consistent
  with Layer 1.

- [x] V-detector (variable-radius) chosen over fixed-radius; rationale recorded.
- [x] PCA-reduced variant included as a measured comparison.
- [ ] Track 1 implements detector generation + coverage stopping rule.

**Cross-track hand-off:** Track 2 sweeps detector count, margin, and PCA dimension.

---

## 3. Self-tolerance calibration — bound the self false-flag rate

**Status:** DECIDED — see P3-D3.

The "calibrated" half. Set the repertoire (detector count / radius margin / decision
threshold) so the **held-out-self false-flag rate ≤ target α** (e.g. 0.05). In immune
terms this *is* self-tolerance — negative selection deletes self-reactive detectors to
keep autoimmunity (false flags on self) low; α is the tolerated self-reactivity. Yields:

- a **calibrated binary** self/non-self decision with a checked error rate (parallels
  Layer-1's EVT p-value, different output), and
- a **continuous non-self score** (e.g. depth of detector coverage / nearest-detector
  margin) for ranking and for the AUROC in deliverable 4.

α is set on a held-out split of self, never on the queries (leakage-safe).

- [x] Self-tolerance = held-out-self false-flag ≤ α; α the immune self-reactivity bound.
- [x] Both a calibrated binary decision and a continuous non-self score emitted.
- [ ] Track 1 implements the calibration on a held-out-self split.

**Cross-track hand-off:** Track 2 sweeps α and reports the coverage/error trade-off curve.

---

## 4. Evaluation + Layer-1 convergence

**Status:** DECIDED — see P3-D4.

Evaluated on the **same benchmark** as Layer 1, so the two are directly comparable:

- **Held-out-family AUROC** (reuses the P2-D4 protocol on the layer-22 reference): withhold
  a Pfam-35 family → its members are non-self queries → AUROC of the non-self score
  separating them from retained self. Distribution over the 903 ≥5-member families.
- **Calibration check:** does the α self-tolerance bound actually hold on held-out self?
- **Convergence with Layer 1 (the load-bearing result):** Spearman rank-correlation between
  the negative-selection non-self score and the Layer-1 EVT novelty on shared queries, plus
  **top-K flag overlap** (Jaccard) — do the two independent methods surface the *same*
  high-novelty genes? Agreement = convergent evidence (the §4 physics "independent lines
  converge" principle), which strengthens the Layer-1 headline.
- **Honest stance:** if negative selection *underperforms* Layer-1's AUROC, that is reported
  plainly. The Layer-3 contribution is the biologically-authentic calibrated self/non-self
  discriminator + the convergence result — **not** beating Layer 1. We are not betting the
  phase on it winning.

- [x] Held-out-family AUROC + calibration + EVT-convergence (Spearman + top-K Jaccard) fixed.
- [x] Underperformance is an acceptable, reportable outcome.
- [ ] Track 1 builds the eval; Track 2 scales/sweeps.

**Cross-track hand-off:** Track 2 produces the comparison tables/figures for review.

---

## 5. Frozen hyperparameters + immune narrative

**Status:** DECIDED (declared pre-result) — see P3-D5.

Declared before any Layer-3 result, not tuned on it (same anti-fishing rule as P2-D5):

- **Detectors:** candidate count (e.g. 10,000) and coverage stopping rule; radius **margin**;
  radius **floor**.
- **Self-tolerance α = 0.05** (headline), sensitivity over {0.01, 0.05, 0.10} as appendix.
- **PCA variant:** dims ∈ {full-1280, 100, 50}, retained variance reported.
- **Convergence:** top-K overlap at K ∈ {50, 100, 500, 1000} (matches P2-D5).
- **Seeds:** all detector sampling seeded and recorded.

**Immune narrative / literature:** negative selection (Forrest et al., 1994), V-detector
(Ji & Dasgupta), innate self/non-self recognition (blueprint §5), and the anomaly-detection
lineage the blueprint ties to the author's prior work. Foregrounds *why* distance-from-self
is a principled novelty signal, not ad hoc.

- [x] All hyperparameters declared pre-result.
- [ ] Track 2 commits `config/immune.yaml` before the first run.

**Cross-track hand-off:** Track 2 commits `config/immune.yaml` with these values before the
first Layer-3 run.

---

**Phase 3 Track 1 — Layer-3 method fixed as calibrated negative selection (P3-D1…P3-D5).**
Track 1 builds `darkmatter/immune/` (V-detector + self-tolerance calibration + scoring +
convergence eval) over the existing layer-22 reference embeddings; Track 2 runs the sweeps
and scales. Next Track-1 step: implement the detector + calibration and run the held-out-family
+ convergence eval, then fold the result into the manuscript narrative (Phase 4).
