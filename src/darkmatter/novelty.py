"""Layer 1 calibrated novelty scorer: EVT (GPD) on kNN cosine distance,
per Track 1's P2-D1/P2-D2 spec.

Score: for a query x, d(x) = mean cosine distance to its k nearest
characterised-at-T0 reference neighbours (mean-of-k, not distance-to-kth
-- robust to one spuriously close reference hit). Calibrated on the
reference's own leave-one-out kNN distances (each reference protein's
kNN distance to its k nearest OTHER reference proteins, so same-family
neighbours are included -- the natural known-to-known spread the tail
model needs) via a Generalized Pareto peaks-over-threshold fit above
threshold u: p(x) = P(a characterised protein sits this far out) = the
checked error rate; novelty s(x) = 1 - p(x).

Raw kNN mean-distance is the uncalibrated baseline the EVT layer must
beat (P2-D1) -- HiFi-NN's hand-set 0.38 cutoff is exactly this family of
method; the contribution is the calibration layer on identical geometry,
not a new representation.

k and the GPD threshold quantile are NOT yet locked by Track 1 (see
docs/problems_and_decisions.md) -- callers must treat results as
provisional until those are fixed, per blueprint D5's no-post-hoc-fishing
rule.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def l2_normalize(embeddings: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(embeddings, axis=1, keepdims=True)
    return embeddings / np.clip(norm, 1e-12, None)


def _topk_mean_distance(query_block: np.ndarray, references: np.ndarray, k: int) -> np.ndarray:
    """Mean cosine distance from each query in the block to its k nearest
    references. Assumes both are L2-normalized (cosine similarity = dot
    product, so distance = 1 - similarity)."""
    sim = query_block @ references.T  # (batch, n_refs)
    k = min(k, sim.shape[1])
    topk_sim = np.partition(sim, -k, axis=1)[:, -k:]
    return (1.0 - topk_sim).mean(axis=1)


def knn_mean_distance(queries: np.ndarray, references: np.ndarray, k: int, batch_size: int = 2000) -> np.ndarray:
    """Mean cosine distance from each query to its k nearest references.
    Batched so the full (n_queries, n_refs) similarity matrix is never
    materialized at once."""
    out = np.empty(len(queries), dtype=np.float64)
    for start in range(0, len(queries), batch_size):
        end = min(start + batch_size, len(queries))
        out[start:end] = _topk_mean_distance(queries[start:end], references, k)
    return out


def leave_one_out_distances(references: np.ndarray, k: int, batch_size: int = 2000) -> np.ndarray:
    """Each reference's mean cosine distance to its k nearest OTHER
    references (excludes itself) -- the calibration distribution the GPD
    tail is fit to."""
    n = len(references)
    out = np.empty(n, dtype=np.float64)
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        block = references[start:end]
        sim = block @ references.T  # (batch, n)
        rows = np.arange(start, end)
        sim[np.arange(end - start), rows] = -np.inf  # exclude self
        k_eff = min(k, n - 1)
        topk_sim = np.partition(sim, -k_eff, axis=1)[:, -k_eff:]
        out[start:end] = (1.0 - topk_sim).mean(axis=1)
    return out


def fit_gpd_tail(calibration_distances: np.ndarray, threshold_quantile: float = 0.9) -> dict:
    """Peaks-over-threshold GPD fit to the upper tail of the leave-one-out
    reference distances."""
    u = float(np.quantile(calibration_distances, threshold_quantile))
    exceedances = calibration_distances[calibration_distances > u] - u
    if len(exceedances) < 10:
        raise ValueError(f"only {len(exceedances)} exceedances above u={u:.4f}, too few for a stable GPD tail fit")
    shape, _loc, scale = stats.genpareto.fit(exceedances, floc=0)
    zeta_u = len(exceedances) / len(calibration_distances)  # P(X > u) empirically
    return {
        "u": u,
        "shape": float(shape),
        "scale": float(scale),
        "zeta_u": float(zeta_u),
        "n_exceedances": int(len(exceedances)),
        "n_calibration": int(len(calibration_distances)),
    }


def evt_novelty_score(distances: np.ndarray, gpd: dict) -> np.ndarray:
    """p(x) = P(a characterised protein sits this far out); novelty
    s(x) = 1 - p(x). Points at or below u never entered the tail model,
    so they get p(x) = 1 (novelty 0) by construction -- unremarkable."""
    exceed = np.clip(distances - gpd["u"], 0, None)
    sf = stats.genpareto.sf(exceed, gpd["shape"], loc=0, scale=gpd["scale"])
    p = np.where(distances > gpd["u"], gpd["zeta_u"] * sf, 1.0)
    return 1.0 - p


def auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """AUROC via the Mann-Whitney U / rank-sum identity -- avoids adding
    scikit-learn as a dependency for one metric."""
    labels = np.asarray(labels, dtype=bool)
    n_pos, n_neg = int(labels.sum()), int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        raise ValueError(f"AUROC undefined with n_pos={n_pos}, n_neg={n_neg}")
    ranks = stats.rankdata(scores)
    sum_ranks_pos = ranks[labels].sum()
    return float((sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def precision_at_k(scores: np.ndarray, labels: np.ndarray, k: int) -> float:
    """Fraction of the top-k ranked-by-score items that are true
    positives. Ties (common here -- most queries score exactly 0, see
    P2 diagnostics) are broken by whatever order argsort returns them in,
    not a hidden secondary key."""
    order = np.argsort(-np.asarray(scores), kind="stable")
    k = min(k, len(order))
    return float(np.asarray(labels)[order[:k]].mean())
