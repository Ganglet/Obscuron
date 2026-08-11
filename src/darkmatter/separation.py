"""Build a small labeled subset from real Pfam-hit proteins and measure
whether an embedding model separates known functional categories — the
Phase 1 deliverable in blueprint section 8 ("run the Genos-m/ESM-2
comparison on a small labeled subset and report the result to Track 1").

"Labeled" here means unambiguous: a protein counts as a member of family
X only if it hits exactly one Pfam-35 family at GA threshold, so the
label itself isn't in question — only whether the embedding space
actually separates it from other families.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np


def select_labeled_subset(
    hits: dict[str, list[str]],
    min_family_size: int = 3,
    per_family: int = 4,
    max_families: int | None = None,
) -> dict[str, list[str]]:
    """{pfam_accession: [sequence_ids]} for families with enough unambiguous single-hit members."""
    by_family: dict[str, list[str]] = defaultdict(list)
    for seq_id, pfams in hits.items():
        if len(pfams) == 1:
            by_family[pfams[0]].append(seq_id)

    families = {fam: ids[:per_family] for fam, ids in by_family.items() if len(ids) >= min_family_size}
    if max_families is not None and len(families) > max_families:
        ranked = sorted(families, key=lambda f: -len(by_family[f]))[:max_families]
        families = {f: families[f] for f in ranked}
    return families


def separation_score(embeddings: np.ndarray, labels: list[str]) -> dict:
    """Mean within-family vs across-family cosine similarity — a bigger gap
    means the embedding space actually separates the labeled categories,
    not just noise."""
    norm = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)
    sim = norm @ norm.T

    labels_arr = np.array(labels)
    same_family = labels_arr[:, None] == labels_arr[None, :]
    not_diagonal = ~np.eye(len(labels), dtype=bool)

    within = sim[same_family & not_diagonal]
    across = sim[~same_family & not_diagonal]

    return {
        "n_sequences": len(labels),
        "n_families": len(set(labels)),
        "within_family_mean_cosine": float(within.mean()),
        "across_family_mean_cosine": float(across.mean()),
        "separation_gap": float(within.mean() - across.mean()),
    }
