"""Calibrated V-detector negative selection over embeddings (Layer 3).

Self = characterised reference. A repertoire of variable-radius detectors covers
the non-self complement: candidate centres whose radius (distance to nearest self,
minus a margin) stays above a floor are kept; self-covering candidates are deleted
(the negative-selection step). A query's non-self score is the deepest coverage by
any surviving detector; the decision threshold is calibrated so the held-out-self
false-flag rate <= alpha (self-tolerance). Distinct from Layer 1: a repertoire
covering the complement + a calibrated binary decision, not a kNN distance rank.

Everything works in cosine geometry on L2-normalised vectors (optionally PCA-reduced
first, then re-normalised), matching P1-D13 / P2-D9. See P3-D1..D5.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _l2(x: np.ndarray) -> np.ndarray:
    return x / np.clip(np.linalg.norm(x, axis=1, keepdims=True), 1e-12, None)


def _nearest_self_dist(points: np.ndarray, self_pts: np.ndarray, chunk: int = 2048) -> np.ndarray:
    """Min cosine distance from each point to any self point (unit vectors)."""
    out = np.empty(len(points))
    for a in range(0, len(points), chunk):
        d = 1.0 - points[a:a + chunk] @ self_pts.T
        out[a:a + chunk] = d.min(axis=1)
    return out


def fit_pca(X: np.ndarray, dims: int) -> tuple[np.ndarray, np.ndarray]:
    """Top-`dims` PCA via the covariance eigendecomposition (fast, reusable across folds)."""
    mean = X.mean(axis=0)
    xc = X - mean
    w, v = np.linalg.eigh(xc.T @ xc)          # ascending eigenvalues
    comp = v[:, ::-1][:, :dims].T             # top-dims components, (dims, d)
    return mean, comp


def _self_nn_scale(S: np.ndarray, chunk: int = 2048) -> float:
    """Median nearest-*other*-self cosine distance — the natural length scale."""
    nn = np.empty(len(S))
    for a in range(0, len(S), chunk):
        d = 1.0 - S[a:a + chunk] @ S.T
        rows = np.arange(d.shape[0])
        d[rows, a + rows] = np.inf  # exclude self
        nn[a:a + chunk] = d.min(axis=1)
    return float(np.median(nn))


@dataclass
class VDetector:
    """Variable-radius negative-selection repertoire with self-tolerance calibration.

    margin/floor are expressed as multiples of the self nearest-neighbour scale, so
    the config is geometry-agnostic (full-space or PCA-reduced).
    """
    n_detectors: int = 5000
    margin_factor: float = 1.0
    floor_factor: float = 0.1
    expand: float = 0.10
    sampling: str = "jitter"     # "jitter" (outward-jittered self points) or "box" (bounding box)
    jitter_scale: float = 3.0    # jitter stddev in units of the self NN scale
    pca_dims: int = 50           # 0 = full embedding space
    max_candidate_mult: int = 30
    seed: int = 42

    # fitted state
    centers_: np.ndarray | None = None
    radii_: np.ndarray | None = None
    pca_mean_: np.ndarray | None = None
    pca_comp_: np.ndarray | None = None
    threshold_: float = 0.0
    n_kept_: int = 0

    def _to_space(self, x: np.ndarray) -> np.ndarray:
        if self.pca_comp_ is not None:
            x = (x - self.pca_mean_) @ self.pca_comp_.T
        return _l2(x)

    def _fit_pca(self, X: np.ndarray) -> None:
        if not self.pca_dims or self.pca_comp_ is not None:
            return  # PCA disabled, or already supplied (reused across held-out folds)
        self.pca_mean_, self.pca_comp_ = fit_pca(X, self.pca_dims)

    def fit(self, self_embeddings: np.ndarray) -> "VDetector":
        rng = np.random.default_rng(self.seed)
        self._fit_pca(self_embeddings)
        s = self._to_space(self_embeddings)

        scale = _self_nn_scale(s)
        margin, floor = self.margin_factor * scale, self.floor_factor * scale
        lo, hi = s.min(axis=0), s.max(axis=0)
        span = hi - lo

        centers, radii = [], []
        attempts, cap = 0, self.n_detectors * self.max_candidate_mult
        while len(centers) < self.n_detectors and attempts < cap:
            if self.sampling == "jitter":  # place candidates in the shell around self
                base = s[rng.integers(0, len(s), size=4096)]
                cand = _l2(base + rng.normal(0.0, self.jitter_scale * scale, size=base.shape))
            else:
                cand = _l2(rng.uniform(lo - self.expand * span, hi + self.expand * span,
                                       size=(4096, s.shape[1])))
            r = _nearest_self_dist(cand, s) - margin
            keep = r > floor
            centers.extend(cand[keep])
            radii.extend(r[keep])
            attempts += 4096

        self.centers_ = np.asarray(centers[: self.n_detectors])
        self.radii_ = np.asarray(radii[: self.n_detectors])
        self.n_kept_ = len(self.centers_)
        return self

    def score(self, queries: np.ndarray, chunk: int = 1024, mode: str = "count") -> np.ndarray:
        """Non-self score. mode='count': fraction of detectors covering q (robust,
        higher = more non-self); mode='depth': deepest relative coverage by any detector."""
        q = self._to_space(queries)
        out = np.zeros(len(q))
        for a in range(0, len(q), chunk):
            d = 1.0 - q[a:a + chunk] @ self.centers_.T
            if mode == "count":
                out[a:a + chunk] = (d < self.radii_[None, :]).mean(axis=1)
            else:
                depth = (self.radii_[None, :] - d) / self.radii_[None, :]
                out[a:a + chunk] = np.clip(depth, 0.0, None).max(axis=1)
        return out

    def calibrate(self, heldout_self: np.ndarray, alpha: float) -> float:
        """Set the threshold so at most alpha of held-out self is flagged non-self."""
        self.threshold_ = float(np.quantile(self.score(heldout_self), 1.0 - alpha))
        return self.threshold_

    def flag(self, queries: np.ndarray) -> np.ndarray:
        return self.score(queries) > self.threshold_
