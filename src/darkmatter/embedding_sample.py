"""Build the Stage 3 embedding sample (P1-D6): all positives, a capped
dark-negative query universe, and a phylum-stratified characterised-at-T0
reference -- the ~50-80k budget any novelty scorer (EVT or density-based)
needs as its raw input, independent of which one Track 1 picks.

Three categories, drawn from panel_protein_labels.csv:
- positive: dark-at-T0 AND characterised-T1-proxy -- never subsampled, rare.
- dark_negative: dark-at-T0 but still dark -- the query universe Precision@K
  ranks over. Capped, since it's most of the panel.
- characterised_at_t0: has a Pfam-35 hit -- the "known boundary" a novelty
  score measures distance from. Capped and phylum-stratified so it isn't
  dominated by whichever phyla happen to have the most genomes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def build_embedding_sample(
    labels: pd.DataFrame,
    panel: pd.DataFrame,
    dark_negative_cap: int = 30000,
    characterised_cap: int = 30000,
    seed: int = 42,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    positives = labels[labels["positive_proxy"]].copy()
    positives["category"] = "positive"

    dark_negative_pool = labels[labels["dark_at_t0"] & ~labels["positive_proxy"]]
    n_dark_neg = min(len(dark_negative_pool), dark_negative_cap)
    dark_negatives = dark_negative_pool.sample(n=n_dark_neg, random_state=rng).copy()
    dark_negatives["category"] = "dark_negative"

    characterised_pool = labels[~labels["dark_at_t0"]].merge(
        panel[["accession", "phylum"]], left_on="genome_accession", right_on="accession", how="left"
    )
    n_phyla = characterised_pool["phylum"].nunique()
    per_phylum_cap = max(1, characterised_cap // max(1, n_phyla))
    characterised = characterised_pool.groupby("phylum", group_keys=False).apply(
        lambda g: g.sample(n=min(len(g), per_phylum_cap), random_state=rng)
    )
    if len(characterised) > characterised_cap:
        characterised = characterised.sample(n=characterised_cap, random_state=rng)
    characterised = characterised.drop(columns=["accession", "phylum"], errors="ignore").copy()
    characterised["category"] = "characterised_at_t0"

    sample = pd.concat([positives, dark_negatives, characterised], ignore_index=True)
    return sample[["genome_accession", "protein_id", "category"]]
