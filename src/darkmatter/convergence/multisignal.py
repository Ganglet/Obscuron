"""Genomic-context novelty + multi-signal convergence for Layer 4.

Genomic-context novelty (the second, embedding-independent axis): a dark gene
surrounded by *characterised* genes sits in a well-understood genomic
neighbourhood, so guilt-by-association makes it more tractable (lower context
novelty); a dark gene embedded among *other dark* genes sits in an unexplored
neighbourhood (higher context novelty). Computed from GTDB gene-order, which is
encoded in the protein IDs (`<contig>_<gene-ordinal>`) -- no new model, no
embedding, so it is a genuinely independent signal from the ESM-2 EVT score.

Convergence follows the SS4 principle: rank by each axis, and treat a gene as a
high-confidence novel candidate only when it is high on *both* -- isolated in
sequence space AND in a dark neighbourhood.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _split_contig_ordinal(protein_id: str) -> tuple[str, int] | None:
    """`CP002988.1_390` -> ("CP002988.1", 390). None if the tail is not an
    integer ordinal (cannot be placed on the contig)."""
    contig, _, tail = protein_id.rpartition("_")
    if not contig or not tail.isdigit():
        return None
    return contig, int(tail)


def genomic_context_novelty(
    labels: pd.DataFrame, window: int, min_neighbors: int
) -> dict[str, float]:
    """protein_id -> context novelty = fraction of its on-contig neighbours that
    are dark-at-T0 (i.e. 1 - characterised-fraction). Uses ALL genes on the
    contig as neighbours (not only scored queries), so the neighbourhood is the
    true genomic one. Genes with fewer than `min_neighbors` defined neighbours
    (contig ends, tiny contigs) are omitted -- context is undefined, not zero.

    `labels` needs columns: protein_id, dark_at_t0 (bool-like).
    """
    dark = labels["dark_at_t0"].astype(str).eq("True").to_numpy()
    parsed = [_split_contig_ordinal(p) for p in labels["protein_id"]]

    # group gene indices by contig, ordered by ordinal
    contigs: dict[str, list[tuple[int, int]]] = {}
    for row_idx, pc in enumerate(parsed):
        if pc is None:
            continue
        contig, ordinal = pc
        contigs.setdefault(contig, []).append((ordinal, row_idx))

    ids = labels["protein_id"].to_numpy()
    out: dict[str, float] = {}
    for genes in contigs.values():
        genes.sort()
        ordinals = [o for o, _ in genes]
        rows = [r for _, r in genes]
        dark_flags = dark[rows].astype(np.int8)
        n = len(genes)
        for i in range(n):
            # neighbours within +/- window in gene-order, excluding self
            lo, hi = max(0, i - window), min(n, i + window + 1)
            neigh = np.concatenate([dark_flags[lo:i], dark_flags[i + 1:hi]])
            if neigh.size < min_neighbors:
                continue
            out[ids[rows[i]]] = float(neigh.mean())  # fraction of neighbours dark
    return out


_AA20 = "ACDEFGHIKLMNPQRSTVWY"


def aa_composition(seq: str) -> np.ndarray:
    """20-dim amino-acid frequency vector (fraction of each standard AA over the
    full length, so a high non-standard-character fraction shows up as a deficit).
    This is the raw material for the third, model-free statistical novelty axis
    (blueprint Layer 5 / the SS4 computational-linguistics parallel)."""
    seq = seq.upper()
    n = max(len(seq), 1)
    return np.array([seq.count(a) for a in _AA20], dtype=float) / n


def mahalanobis_novelty(query_comp: np.ndarray, ref_comp: np.ndarray, ridge: float = 1e-6) -> np.ndarray:
    """Compositional novelty = Mahalanobis distance of each query's AA composition
    from the characterised reference's composition distribution (mean + covariance,
    ridge-stabilised). Higher = more compositionally unusual = more novel. Shares no
    input with either the ESM-2 embedding or the genomic context -> a third
    independent axis."""
    mu = ref_comp.mean(0)
    cov = np.cov(ref_comp, rowvar=False) + ridge * np.eye(ref_comp.shape[1])
    inv = np.linalg.inv(cov)
    diff = query_comp - mu
    return np.sqrt(np.clip(np.einsum("ij,jk,ik->i", diff, inv, diff), 0, None))


def high_confidence_convergent_multi(axes: list[np.ndarray], quantile: float) -> np.ndarray:
    """Boolean mask: above the per-axis `quantile` cut on ALL axes (the N-way
    generalisation of high_confidence_convergent)."""
    mask = np.ones(len(axes[0]), dtype=bool)
    for a in axes:
        mask &= a >= np.quantile(a, quantile)
    return mask


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    ra, rb = ra - ra.mean(), rb - rb.mean()
    denom = np.sqrt((ra**2).sum() * (rb**2).sum())
    return float((ra * rb).sum() / denom) if denom > 0 else float("nan")


def _top_jaccard(a: np.ndarray, b: np.ndarray, k: int) -> float:
    """Jaccard overlap of the two top-k sets (highest values)."""
    k = min(k, len(a))
    ta = set(np.argsort(a)[::-1][:k])
    tb = set(np.argsort(b)[::-1][:k])
    return len(ta & tb) / len(ta | tb)


def convergence_stats(evt: np.ndarray, context: np.ndarray, top_ks: list[int]) -> dict:
    """Independence + agreement of the two axes on the shared query set."""
    return {
        "n": int(len(evt)),
        "spearman": _spearman(evt, context),
        "top_jaccard": {k: _top_jaccard(evt, context, k) for k in top_ks},
    }


def high_confidence_convergent(
    evt: np.ndarray, context: np.ndarray, quantile: float
) -> np.ndarray:
    """Boolean mask: high on BOTH axes (>= the per-axis `quantile` cut). This is
    the SS4 convergence set -- confident novelty only where both agree."""
    return (evt >= np.quantile(evt, quantile)) & (context >= np.quantile(context, quantile))
