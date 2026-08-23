"""kNN cosine distance -- brute force, batched to avoid materializing the
full (n_query, n_reference) similarity matrix at once. Interface fixed by
Track 1's docs/Track2_Phase2_scoring_handoff.md.
"""

from __future__ import annotations

import numpy as np


def l2_normalize(embeddings: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(embeddings, axis=1, keepdims=True)
    return embeddings / np.clip(norm, 1e-12, None)


def knn_distance(
    query: np.ndarray,
    reference: np.ndarray,
    k: int = 5,
    reduce: str = "mean",
    exclude_self: bool = False,
    batch_size: int = 2000,
) -> np.ndarray:
    """Assumes query and reference are already L2-normalised (cosine
    similarity = dot product, distance = 1 - similarity).

    exclude_self=True is the leave-one-out calibration pass, where query
    IS reference: each row's own zero-distance self-match is dropped
    before taking the k nearest, so same-family neighbours still count
    (the natural known-to-known spread the tail model needs) but a
    protein never neighbours itself.
    """
    if reduce != "mean":
        raise NotImplementedError(f"reduce={reduce!r} not implemented -- P2-D2 fixes mean-of-k")

    n = len(query)
    out = np.empty(n, dtype=np.float64)
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        block = query[start:end]
        sim = block @ reference.T  # (batch, n_ref)

        if exclude_self:
            rows = np.arange(start, end)
            sim[np.arange(end - start), rows] = -np.inf

        k_eff = min(k, sim.shape[1] - (1 if exclude_self else 0))
        topk_sim = np.partition(sim, -k_eff, axis=1)[:, -k_eff:]
        out[start:end] = (1.0 - topk_sim).mean(axis=1)
    return out
