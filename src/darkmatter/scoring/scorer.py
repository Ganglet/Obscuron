"""NoveltyScorer -- ties knn.py + gpd.py together per the interface fixed
in docs/Track2_Phase2_scoring_handoff.md. Identical across embedding arms
(P2-D2); only the input matrix changes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from darkmatter.scoring.gpd import GPDTail, fit_gpd_tail, tail_pvalue
from darkmatter.scoring.knn import knn_distance


class NoveltyScorer:
    def __init__(self, config: dict):
        scorer_cfg = config["scorer"]
        self.metric = scorer_cfg["metric"]
        if self.metric != "cosine":
            raise NotImplementedError(f"metric={self.metric!r} not implemented -- only cosine (P2-D2)")
        self.k = scorer_cfg["knn"]["k"]
        self.reduce = scorer_cfg["knn"]["reduce"]
        self.threshold_quantile = scorer_cfg["gpd"]["threshold_quantile"]

        self._reference: np.ndarray | None = None
        self._loo_distances: np.ndarray | None = None
        self._tail: GPDTail | None = None

    def fit(self, reference: np.ndarray) -> None:
        self._reference = reference
        self._loo_distances = knn_distance(reference, reference, k=self.k, reduce=self.reduce, exclude_self=True)
        self._tail = fit_gpd_tail(self._loo_distances, threshold_quantile=self.threshold_quantile)

    @property
    def tail(self) -> GPDTail:
        if self._tail is None:
            raise RuntimeError("call fit() before accessing tail")
        return self._tail

    @property
    def loo_distances(self) -> np.ndarray:
        if self._loo_distances is None:
            raise RuntimeError("call fit() before accessing loo_distances")
        return self._loo_distances

    def score(self, query: np.ndarray) -> pd.DataFrame:
        if self._tail is None or self._reference is None:
            raise RuntimeError("call fit() before score()")
        distance = knn_distance(query, self._reference, k=self.k, reduce=self.reduce, exclude_self=False)
        pvalue = tail_pvalue(self._tail, distance, ref_distances=self._loo_distances)
        return pd.DataFrame({"distance": distance, "pvalue": pvalue, "novelty": 1.0 - pvalue})

    def score_baseline(self, query: np.ndarray) -> np.ndarray:
        """Raw kNN distance only, no GPD -- the HiFi-NN-style baseline the
        calibrated headline must beat (P2-D1)."""
        if self._reference is None:
            raise RuntimeError("call fit() before score_baseline()")
        return knn_distance(query, self._reference, k=self.k, reduce=self.reduce, exclude_self=False)
