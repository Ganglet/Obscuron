# Phase 2 Track 1 — Calibrated Novelty Scorer Design

**Phase:** 2 — Calibrated Novelty Scoring (Layer 1)
**Track:** 1 — Methodology, Design & Analysis
**Status:** In progress
**Branch:** `phase-2-track-1`

---

## Objective

Turn the Phase-1 benchmark and embeddings into the project's core method: a
*calibrated* novelty score for dark-at-T₀ genes, with an error rate that is
**checked rather than asserted**, validated on the leakage-controlled
retrospective benchmark. Track 1 fixes the scorer, the reference, the grading,
and the frozen hyperparameters; Track 2 implements them against this spec.
Nothing about the ranking is tuned on results. The headline the title commits
to — *calibrated* novelty — is delivered here or nowhere.

Five deliverables, in dependency order:

1. Scorer family — the calibrated headline vs the baseline it must beat.
2. Reference set + score construction — what "characterised space" is and how distance becomes the raw score.
3. Evaluation wiring — how the retrospective positives and the held-out-family test feed the three frozen metrics.
4. Held-out-family protocol — the controlled novelty ground truth for AUROC + calibration.
5. Frozen hyperparameters — declared before any result, anti-fishing.

Blueprint references: §6 Layer 1, §7 (recommended approach, metrics, go/no-go).
Decisions: **P2-D1–P2-D5**, building on **P1-D3** (labels), **P1-D7** (leakage),
**P1-D11** (novelty claim), **P1-D13** (embedding geometry).

---

## 1. Scorer family

**Status:** DECIDED — see P2-D1.

**Headline = EVT peaks-over-threshold (GPD) on the kNN cosine-distance score.
Baseline = raw kNN distance (the HiFi-NN-style heuristic cutoff, P1-D11).**

The baseline *ranks* dark genes but reports no honest error rate — it is exactly
the family of method the closest prior art (HiFi-NN's hand-set 0.38 cutoff) uses.
The headline fits a Generalized Pareto Distribution to the **tail** of the
known-to-known distance distribution (Pickands–Balkema–de Haan: threshold
exceedances converge to a GPD), turning the *same* distance into a **calibrated
tail p-value with a checked error rate**. The contribution is precisely the
calibration layer on identical geometry — not a new representation.

Per-class EVM / OpenMax (a Weibull tail per Pfam family) was **rejected as the
headline**: the panel has thousands of families, many with a handful of members
→ thin, unstable per-class tails. It is retained as a **secondary comparison**
only if family counts allow (§5). We score distance from *all* known space, which
is what the single global GPD measures data-efficiently.

**Score definition.**
- `d(x)` = mean cosine distance of query `x` to its `k` nearest
  characterised-reference neighbours (mean-of-`k`; `k` frozen in §5).
- Calibrate on the reference's own **leave-one-out** kNN distances `{d(r)}`; fit a
  GPD to the tail above threshold `u`.
- For `d(x) > u`: tail p-value
  `p(x) = ζ_u · (1 − G(d(x) − u; ξ, β))`, where `ζ_u` = fraction of the reference
  exceeding `u`. `p(x)` = estimated **P(a characterised protein sits this far
  out)** = the calibrated false-flag rate at that cutoff.
- Novelty `s(x) = 1 − p(x)`; for `d(x) ≤ u`, `s` from the empirical CDF (the tail
  is the region of interest — that is where "novel" lives).

- [x] Headline vs baseline fixed (P2-D1).
- [x] Score definition (kNN → GPD tail p-value) fixed.
- [x] EVM/OpenMax scoped as secondary, gated on family counts.
- [ ] Track 2 implements + first scored run (ESM-2 arm).

**Cross-track hand-off:** Track 2 builds a `darkmatter/scoring/` module — kNN over
the L2-normalised embeddings (P1-D13), leave-one-out reference calibration, GPD
fit (`scipy.stats.genpareto` MLE), returning `p(x)` and `s(x)`. Baseline = the
same kNN distance without the GPD step. EVM/OpenMax secondary, gated on §5.

---

## 2. Reference set + score construction

**Status:** DECIDED — see P2-D2.

- **Reference (= "characterised space") = characterised-at-T₀ panel proteins**
  (a Pfam-35 GA-threshold hit, P1-D3). **Queries = dark-at-T₀ panel proteins**
  (no Pfam-35 hit) — the population that gets scored. The 4,138 **positives** are
  the subset of queries later characterised by T₁ (P1-D12).
- **Geometry = cosine** (L2-normalised embeddings, P1-D13); distance
  `= 1 − cosine similarity`. No re-embedding; the scorer inherits Phase-1 vectors.
- **kNN score = mean of the `k` nearest characterised distances** (mean-of-`k`,
  not distance-to-`k`th alone — robust to a single spuriously close reference hit;
  standard for kNN novelty).
- **Reference calibration = leave-one-out:** each reference protein's kNN distance
  excludes itself but includes same-family neighbours — the natural known-to-known
  spread the tail model must reflect.
- **Per-arm leakage (P1-D7):** the scorer is **identical across arms; only the
  embedding changes.** ESM-2 = leakage-clean headline; Genos-m = leakage-controlled
  (positives restricted post-R220, plus the full-window number), reported
  separately, never pooled.

- [x] Reference = characterised-at-T₀; queries = dark-at-T₀ (P1-D3).
- [x] Cosine geometry inherited (P1-D13); mean-of-`k` score fixed.
- [x] Leave-one-out reference calibration fixed.
- [ ] Track 2 builds reference/query matrices from Stage-3 embeddings + labels.

**Cross-track hand-off:** Track 2 assembles the reference and query embedding
matrices from the Stage-3 output + the dark/positive labels; `k` and the distance
metric are config, not hard-coded.

---

## 3. Evaluation wiring

**Status:** DECIDED — see P2-D3.

The three frozen metrics (blueprint §7), each fed by a **specific** ground truth —
and calibration is certified on the *controlled* ground truth, **never** the
selection-biased positives.

- **Precision@K (ecological validity).** Rank dark queries by `s(x)`; top-`K`
  precision = fraction that are positives. Report at frozen `K` (§5) as a curve,
  and as **lift over the scored set's own positive base rate** (positives /
  dark queries ranked — ≈1.4% at natural dark-population density, *enriched* in the
  capped embedded evaluation set since P1-D6 caps dark negatives; report the exact
  set value). Lift, not raw precision, is the signal statistic. **Selection-bias
  caveat (P1-D3) stated, not buried:** positives skew *moderate*-novelty (near-known
  genes get characterised first; maximal-novelty genes are too far out to have been
  characterised yet), so Precision@K measures **prioritisation value**, not
  "novelty = characterisability." A monotone lift over base rate is the claim; a
  perfect top-`K` is neither expected nor the target.
- **AUROC (sensitivity).** From the held-out-family construction (§4): separate
  held-out-family members (label **novel**) from retained knowns (label **known**).
  Controlled ground truth, no selection bias. Report the distribution over held-out
  families.
- **Calibration check (honesty — the headline).** Reliability diagram on the
  **held-out-family** ground truth: bin queries by predicted p-value, plot the
  observed novel-fraction per bin against the diagonal. This directly tests whether
  the GPD's "checked error rate" is honest. **Built on held-out-family (controlled),
  not on characterised-by-T₁ (biased)** — the load-bearing methodological choice:
  we never ask the selection-biased positive set to certify calibration.

One-line split: **held-out-family certifies the score is *sensitive and honest*;
the retrospective positives show it surfaces *real, later-confirmed biology*.**

- [x] Precision@K driven by positives, reported as lift over the set base rate.
- [x] AUROC + calibration driven by held-out-family (controlled), not positives.
- [x] Selection-bias framing fixed (prioritisation, not determination).
- [ ] Track 2 implements `evaluate_scorer.py` reading a frozen config (§5).

**Cross-track hand-off:** Track 2 implements `evaluate_scorer.py` — Precision@K +
lift, AUROC over held-out families, reliability diagram — all metrics, `K`, and
seeds read from the frozen `config/scorer.yaml` written **before** results.

---

## 4. Held-out-family protocol

**Status:** DECIDED — see P2-D4.

The controlled novelty ground truth that §3's AUROC and calibration depend on.

- **Construction:** withhold one Pfam-35 family `F` → remove all `F` members from
  the reference → `F`'s members become queries labelled **novel** (with `F` gone
  they should sit far from the remaining known space). Negatives (label **known**)
  = a matched sample of characterised proteins whose family is *not* held out.
- **Matched controls:** match novel vs known on **sequence length, GTDB phylum,
  and family size**, so AUROC reflects novelty detection rather than a
  length/rarity confound.
- **Leakage-safe by construction:** the embedding may still place `F`'s members
  near where `F` sat (it saw `F` in pretraining) — that is exactly the stress test:
  can distance-to-*remaining*-reference flag them anyway? The test isolates the
  **scorer's** sensitivity with the embedding held fixed; it is not a claim about
  the embedding's own leakage (that is P1-D7's job).
- **Sampling:** hold out all Pfam-35 families with ≥ `m` members in the panel
  reference; if that set is large, a phylum/size-stratified sample of `N` families
  (§5), fixed seed. Aggregate → AUROC distribution + pooled calibration curve.

- [x] Withhold → score → aggregate construction fixed.
- [x] Matched controls (length, phylum, family size) fixed.
- [x] Leakage-safety argument recorded (scorer sensitivity, embedding fixed).
- [ ] Track 2 implements the withhold-score-aggregate loop.

**Cross-track hand-off:** Track 2 implements `heldout_family_eval.py` — the
withhold/score/aggregate loop, controls sampled with the frozen seed.

---

## 5. Frozen hyperparameters

**Status:** DECIDED (declared pre-result) — see P2-D5.

Fixed **now**, before any scorer output, and not tuned on results — the
anti-fishing commitment the metric-freeze already made. These defaults are
Track-1's call; change them only by a reviewed edit that states a reason (mirrors
the `config/snapshots.yaml` rule).

- **`k` (kNN):** headline **`k = 5`**; sensitivity reported over `{1, 5, 10, 20}`
  as an appendix, **not** used to pick the headline.
- **GPD tail threshold `u`:** **90th percentile** of the reference leave-one-out
  kNN distances (upper 10% = the tail), with a mean-residual-life /
  parameter-stability plot as a robustness check — the *rule* is fixed, not a
  hand-tuned value.
- **GPD fit:** MLE (`scipy.stats.genpareto`); report shape `ξ`, scale `β`, and a
  tail goodness-of-fit (QQ / Anderson–Darling).
- **Precision@K:** `K ∈ {50, 100, 500, 1000}` + the full Precision-vs-`K` curve;
  lift = Precision@K ÷ (positives ÷ dark queries ranked).
- **Held-out families:** all Pfam-35 families with ≥ **`m = 5`** members in the
  panel reference; if > 500 such families, a phylum/size-stratified sample of 500;
  seed = the Phase-1 panel seed.
- **Fair ESM-2-vs-Genos-m comparison:** restrict to genes ≤ 1024 aa (P1-D13) so
  the length asymmetry does not confound the arm comparison.
- **Seeds:** every sampling step seeded and recorded in the frozen config.

- [x] All values above declared pre-result.
- [ ] Track 2 commits `config/scorer.yaml` with these values before the first run.

**Cross-track hand-off:** Track 2 commits `config/scorer.yaml` with every value
above **before** the first scoring run; the evaluation scripts read only from it.

---

**Phase 2 Track 1 — scorer design fixed (P2-D1–P2-D5).** Track 2 implements
`darkmatter/scoring/` + `evaluate_scorer.py` + `heldout_family_eval.py` against
this spec, ESM-2 arm first (leakage-clean headline). Next Track-1 step: review the
first scored run + calibration curve, then Phase 3 (immune-inspired self/non-self
narrative, reusing these embeddings).
