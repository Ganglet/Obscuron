# Layer 3 Sweeps, Full-Scale Dark-Query Flagging & Detector-Count Saturation

**Phase:** 3 — Immune-Inspired Self/Non-Self Discrimination
**Owner:** Track 2 (Rayyan)
**Branch:** `phase-3-track-2-self-nonself-model`
**Status:** Complete. Closes the Track-2 hand-off from `docs/Track1_phase3_immune_design.md` (detector-count/margin/PCA sweeps, alpha sweep, full-scale dark-query flagging) and P3-D6's open follow-up (negative selection was only ever tested against the reference, never the real dark queries). See P3-D7 and P3-D8 in `docs/problems_and_decisions.md` for the full decision-log entries this doc summarizes and expands on.

## Objective

Track 1 (the user) built and reserved Layer 3's core module (`darkmatter/immune/`, V-detector negative selection with self-tolerance calibration) and delivered its first result (P3-D6: calibration works, detection is a weaker-but-honest complement to Layer 1). Per the design doc's own division of labour, Track 2's half was to sweep the frozen hyperparameters for headroom, scale flagging to the real dark-query population (not just the held-out-family reference test), and report whatever the numbers actually say — including anything that disagrees with Layer 1.

## What Was Done

- **Independent reference re-embedding**: re-embedded the characterised reference (15,360 proteins, layer 22, fp32) on the RTX 4060 rather than reusing P3-D6's `.npz` file, since embedding matrices are gitignored per the project's data-provenance standard (P1-D4) and don't travel across machines. The resulting family structure (903 Pfam-35 families with ≥5 members) matched P3-D6's number exactly — a clean cross-machine reproduction check before trusting anything built on top of it.
- **Dark-query embedding, new**: P3-D6 never had layer-22 embeddings for the dark-query set — only the reference. Embedded the dark queries at the same layer/settings so Layer 3 could finally be scored against real, unlabelled-at-T0 sequences instead of only the held-out-family reference subset.
- **Compute-time scope decision, documented rather than hidden**: the full dark-query population is 34,138 sequences. Embedding all of it plus the reference at fp32 on a single 8GB laptop GPU was projected at ~8.5 hours. Scoped down to all 4,138 positives (never subsampled) plus a stratified 10,000-of-30,000 `dark_negative` sample (seed=42) — 14,138 queries, ~5 hours total. Same documented-tradeoff pattern as P1-D6/P2-D10: a real, stated scope reduction, not a silent claim of full-scale coverage.
- **`scripts/immune_sweep.py`**: checkpointed, resumable sweep driver (same resilience pattern as P1-D10) covering three frozen-config axes on 30 sampled held-out families:
  - **Detector count** `[1000, 5000, 10000, 20000]` — AUROC climbs 0.59 → 0.74 (the frozen production default) → 0.83 → 0.86, still rising at the top tested value. Reported as unexploited headroom, not chased into the frozen config (P2-D5/P3-D5's anti-fishing rule — the config was frozen *before* this result existed).
  - **Margin factor** `[0.5, 1.0, 1.5, 2.0]` — flat (0.627–0.629), not a load-bearing hyperparameter in this range.
  - **PCA dimension** `[0, 50, 100]` — validates P3-D2's original design call: no PCA collapses to near-chance (0.53), PCA-50 (the frozen default) reaches 0.629, PCA-100 is slightly worse (0.617).
- **Alpha sweep** (full 5,000-detector repertoire, real held-out self + real dark queries): self-tolerance holds cleanly at every tested α — actual held-out-self flag rate always well under the nominal target (0.010/0.023/0.060 at α=0.01/0.05/0.10). The α=0.05 point (0.024 vs P3-D6's original 0.024) is a second, independent confirmation the calibration mechanism is sound, not a coincidence of the first run.
- **Found and reported an inversion, not smoothed over**: at α=0.05, Layer 3's flag rate is *enriched* in eventual positives over still-dark queries (20.1% vs 8.6%, a 2.3× ratio) — the opposite direction from Layer 1's Precision@K finding (novelty anti-predicts near-term characterisation). Reported as an honest, unresolved divergence for Track 1 to interpret rather than reconciled or hidden. Plausible mechanism written up in P3-D7: Layer 3's coverage score responds to local geometric gaps near the self manifold, which may place moderately-novel, near-known genes (exactly the population that gets characterised) in detector-coverable pockets that far-dark queries reached by sparser candidate sampling are not.
- **Full-scale Layer 1 / Layer 3 convergence** (all 14,138 real dark queries, not the ~700-member held-out-family test): Spearman ρ = 0.485, clearly stronger than P3-D6's held-out-family-only +0.32. Flagged the accompanying caveat that top-K Jaccard staying low (0.00 at K=50 → 0.12 at K=1000) despite this positive correlation is an artefact of Layer 3's score being heavily tied at zero for most queries (a discrete coverage-count, not a continuous distance) — the Spearman number is the trustworthy statistic here, the Jaccard number is not a like-for-like comparison to Layer 1's continuous ranking.
- **Resilience infrastructure**: `scripts/run_immune_pipeline_resilient.ps1`, a Startup-folder auto-resume launcher — checkpointed steps, `HF_HUB_OFFLINE=1` (both models already cached locally), and a watchdog restart on stall. Paid off directly: the run survived one real multi-hour idle gap (machine idle ~14:56→17:37) and resumed with under a checkpoint's worth of lost work, unattended.
- **P3-D8, follow-up saturation sweep**: the detector-count sweep above was still climbing at its top tested value (20,000 → 0.86), leaving the true ceiling unmeasured. `scripts/immune_detector_saturation.py` extended the same protocol to 30,000 / 40,000 / 50,000 detectors. The 20,000-detector point reproduced exactly (0.862368882894188 both runs — a clean determinism check) before trusting the new points. Marginal AUROC gain collapses sharply past 20-30k: 10k→20k +0.030, 20k→30k +0.027, 30k→40k +0.001, 40k→50k +0.009 (noise-level on a 30-family sample) — a real plateau, not still-rising noise. **Saturation point ≈ 30,000 detectors, AUROC ≈ 0.89**, above the frozen production default's 0.74 but still short of Layer 1's full-scale 0.962.

## Why (Key Decisions)

**Why re-embed the reference instead of reusing P3-D6's file?**
Embedding matrices are gitignored (large binaries, P1-D4) and don't travel across the two machines this project runs on. Re-embedding and getting an exact family-count match first was cheap insurance against silently building the rest of this work on a subtly different reference.

**Why a stratified 14,138-query subset instead of the full 34,138?**
Full-scale fp32 embedding of both populations on one 8GB laptop GPU was projected at ~8.5 hours. All positives were kept in full (the label that actually matters for Precision-style reads); the much larger dark-negative pool was stratified-sampled instead of silently truncated, and the reduction is stated here rather than claimed as full coverage.

**Why not chase the detector-count headroom into the frozen config?**
`config/immune.yaml`'s hyperparameters were frozen *before* any result existed, specifically so a finding like "5,000 undershoots 30,000" can't be chased post-hoc (P2-D5/P3-D5's anti-fishing rule). The measurement is reported so Track 1 can make a reviewed call with real data in hand (P3-D8) — re-freezing the default is explicitly left open as Track 1's decision, not made unilaterally here.

**Why report the positive-enrichment inversion instead of tuning it away?**
It's a genuine, reproducible finding (confirmed at three α values) that disagrees with Layer 1's own selection-bias story. A method that reports the axis that disagrees with its own headline is more trustworthy than one that only reports the axes that agree — the same principle P4-D8 later applied to Genos-m's convergence-axis inversion.

## Outputs

| Output | Description |
|---|---|
| `scripts/reembed_reference_immune.py`, `scripts/reembed_dark_queries.py` | Cross-machine reference and dark-query re-embedding at layer 22, committed. |
| `scripts/immune_sweep.py` | Checkpointed, resumable sweep driver (detector count, margin, PCA, alpha), committed. |
| `scripts/immune_detector_saturation.py` | Follow-up detector-count sweep extending to 50,000, committed. |
| `results/immune_sweep_{detector_count,margin,pca_dims,alpha}.csv` | Sweep results, committed. |
| `results/immune_sweep_detector_saturation.csv` | Saturation-point sweep results (P3-D8), committed. |
| `results/immune_scale_summary.json` | Full-scale flag-rate and convergence summary. |
| `results/immune_dark_query_flags.csv` | Per-query Layer-1 novelty + Layer-3 score/flag, all 14,138 scored queries. |
| `results/figures/immune_sweep_summary.png` | Sweep summary figure. |
| `scripts/run_immune_pipeline_resilient.ps1` | Checkpointed, watchdog-restarted, Startup-folder auto-resume launcher. |

**Next:** reported to Track 1 for interpretation against P3-D6's acceptance sanity (the detector-count headroom and the positive-vs-dark_negative flag-rate divergence from Layer 1); folded into the manuscript's Layer 3 section alongside P3-D6 (Phase 4). Track 1's decision on whether to re-freeze the production detector count at ~30,000 remains open — see P3-D8 and `docs/problems_and_decisions.md`.
