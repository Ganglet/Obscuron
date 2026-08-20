"""EVT peaks-over-threshold GPD tail model (P2-D1). Fits a Generalized
Pareto Distribution to the upper tail of the reference's leave-one-out
kNN distances, turning a raw distance into a calibrated tail p-value:
p(x) = P(a characterised protein sits this far out).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass
class GPDTail:
    u: float
    zeta_u: float
    xi: float
    beta: float
    n_exceedances: int
    n_calibration: int


def fit_gpd_tail(ref_distances: np.ndarray, threshold_quantile: float = 0.90) -> GPDTail:
    u = float(np.quantile(ref_distances, threshold_quantile))
    exceedances = ref_distances[ref_distances > u] - u
    if len(exceedances) < 10:
        raise ValueError(f"only {len(exceedances)} exceedances above u={u:.4f}, too few for a stable GPD tail fit")
    xi, _loc, beta = stats.genpareto.fit(exceedances, floc=0)
    zeta_u = len(exceedances) / len(ref_distances)
    return GPDTail(
        u=u, zeta_u=float(zeta_u), xi=float(xi), beta=float(beta),
        n_exceedances=len(exceedances), n_calibration=len(ref_distances),
    )


def tail_pvalue(model: GPDTail, distances: np.ndarray, ref_distances: np.ndarray | None = None) -> np.ndarray:
    """d > u: calibrated GPD tail p-value -- the headline. d <= u: empirical
    survival fraction from the reference distribution -- not the headline
    claim (the tail is the region of interest), but still a defined value
    rather than a hard-coded 1.0 when a reference distribution is supplied."""
    distances = np.asarray(distances, dtype=np.float64)
    p = np.empty_like(distances)

    above = distances > model.u
    exceed = np.clip(distances[above] - model.u, 0, None)
    sf = stats.genpareto.sf(exceed, model.xi, loc=0, scale=model.beta)
    p[above] = model.zeta_u * sf

    if ref_distances is not None and (~above).any():
        ref_sorted = np.sort(ref_distances)
        ranks = np.searchsorted(ref_sorted, distances[~above], side="left")
        p[~above] = 1.0 - ranks / len(ref_sorted)
    else:
        p[~above] = 1.0

    return p


def novelty(model: GPDTail, distances: np.ndarray, ref_distances: np.ndarray | None = None) -> np.ndarray:
    return 1.0 - tail_pvalue(model, distances, ref_distances)


def diagnostics(ref_distances: np.ndarray, model: GPDTail) -> dict:
    """Mean-residual-life points and QQ vs the fitted GPD (P2-D5). Uses a
    KS test against the fitted GPD as the tail goodness-of-fit stat --
    scipy's Anderson-Darling (stats.anderson) only supports a fixed list
    of named distributions and doesn't take arbitrary genpareto
    parameters, so KS is the available substitute for the same purpose
    (a formal p-value on whether the exceedances plausibly came from the
    fitted tail)."""
    exceedances = ref_distances[ref_distances > model.u] - model.u

    thresholds = np.quantile(ref_distances, np.linspace(0.5, 0.99, 20))
    mrl = []
    for t in thresholds:
        excess = ref_distances[ref_distances > t] - t
        if len(excess) > 1:
            mrl.append({"threshold": float(t), "mean_excess": float(excess.mean()), "n": int(len(excess))})

    n = len(exceedances)
    probs = (np.arange(1, n + 1) - 0.5) / n
    theoretical_q = stats.genpareto.ppf(probs, model.xi, loc=0, scale=model.beta)
    empirical_q = np.sort(exceedances)

    ks_stat, ks_pvalue = stats.kstest(exceedances, "genpareto", args=(model.xi, 0, model.beta))

    return {
        "mean_residual_life": mrl,
        "qq_theoretical": theoretical_q.tolist(),
        "qq_empirical": empirical_q.tolist(),
        "ks_statistic": float(ks_stat),
        "ks_pvalue": float(ks_pvalue),
        "n_exceedances": n,
    }
