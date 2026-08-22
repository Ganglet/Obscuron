"""Phase 2 scorer figures (Track 1 analysis).

Consumes the frozen result schema Track 2 emits
(docs/Track2_Phase2_scoring_handoff.md -> "Result schema"):

  scorer_{arm}_ranked.csv        protein_id, novelty, pvalue, distance, is_positive, phylum, length_aa
  scorer_{arm}_gpd_qq.csv        empirical, theoretical
  heldout_{arm}_auroc.csv        family_id, n_members, phylum, auroc
  heldout_{arm}_calibration.csv  predicted, observed, n

Every figure is a claim in the paper, so each function documents the shape a
*healthy* result takes (the acceptance-sanity in the hand-off doc).
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

HEADLINE = "#3B6FB0"   # GPD novelty
BASELINE = "#C0603A"   # raw kNN distance
REF = "#8A8A8A"        # diagonal / chance


def _ax(ax):
    if ax is None:
        _, ax = plt.subplots(figsize=(4.6, 3.6), constrained_layout=True)
    return ax


def precision_at_k(ranked: pd.DataFrame, score_col: str, k_list) -> tuple[pd.DataFrame, float]:
    """P@K and lift over the query base rate for a ranking by ``score_col`` (descending)."""
    order = ranked.sort_values(score_col, ascending=False)
    pos = order["is_positive"].to_numpy().astype(float)
    base = pos.mean()
    rows = []
    for k in k_list:
        kk = min(int(k), len(pos))
        p = pos[:kk].mean()
        rows.append({"K": kk, "precision": p, "lift": (p / base) if base > 0 else np.nan})
    return pd.DataFrame(rows), base


def precision_at_k_curve(ranked: pd.DataFrame, k_list, ax=None):
    """Headline (GPD novelty) vs baseline (raw kNN distance) lift.

    Healthy: headline lift > 1 and >= baseline at every K (calibration reranks, it
    never hurts). Interpret against the enriched set base rate, not the natural ~1.4%.
    """
    ax = _ax(ax)
    head, base = precision_at_k(ranked, "novelty", k_list)
    basel, _ = precision_at_k(ranked, "distance", k_list)
    ax.plot(head["K"], head["lift"], "o-", color=HEADLINE, label="GPD novelty (headline)")
    ax.plot(basel["K"], basel["lift"], "s--", color=BASELINE, label="raw kNN distance (baseline)")
    ax.axhline(1.0, color=REF, lw=1, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("K")
    ax.set_ylabel("lift over base rate")
    ax.set_title(f"Precision@K lift  (base rate {base:.3f})")
    ax.legend(frameon=False, fontsize=8)
    return ax.figure


def reliability_diagram(calib: pd.DataFrame, ax=None):
    """Calibration on the held-out-family ground truth (the headline honesty claim).

    Healthy: points hug the diagonal -> the GPD p-value's checked error rate is honest.
    """
    ax = _ax(ax)
    ax.plot([0, 1], [0, 1], color=REF, ls=":", lw=1, label="perfect")
    ax.plot(calib["predicted"], calib["observed"], "o-", color=HEADLINE, label="observed")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("predicted novelty")
    ax.set_ylabel("observed novel fraction")
    ax.set_title("Calibration (held-out family)")
    ax.legend(frameon=False, fontsize=8)
    return ax.figure


def auroc_distribution(auroc: pd.DataFrame, ax=None):
    """Sensitivity: per-held-out-family AUROC.

    Healthy: median clearly > 0.5 (the separation gap was 0.036 in Phase 1 -> signal,
    not a landslide).
    """
    ax = _ax(ax)
    vals = auroc["auroc"].to_numpy()
    med = float(np.median(vals))
    ax.hist(vals, bins=20, range=(0.4, 1.0), color=HEADLINE, alpha=0.85)
    ax.axvline(0.5, color=REF, ls=":", lw=1, label="chance")
    ax.axvline(med, color=BASELINE, lw=1.5, label=f"median {med:.3f}")
    ax.set_xlabel("held-out-family AUROC")
    ax.set_ylabel("families")
    ax.set_title("Sensitivity (held-out-family AUROC)")
    ax.legend(frameon=False, fontsize=8)
    return ax.figure


def gpd_qq_plot(qq: pd.DataFrame, ax=None):
    """GPD tail goodness-of-fit. Healthy: points on the diagonal."""
    ax = _ax(ax)
    lo = float(min(qq["empirical"].min(), qq["theoretical"].min()))
    hi = float(max(qq["empirical"].max(), qq["theoretical"].max()))
    ax.plot([lo, hi], [lo, hi], color=REF, ls=":", lw=1)
    ax.plot(qq["theoretical"], qq["empirical"], "o", ms=4, color=HEADLINE)
    ax.set_xlabel("theoretical GPD quantile")
    ax.set_ylabel("empirical quantile")
    ax.set_title("GPD tail fit (QQ)")
    return ax.figure
