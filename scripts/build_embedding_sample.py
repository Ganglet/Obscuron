#!/usr/bin/env python
"""CLI: build the Stage 3 embedding sample manifest (P1-D6) -- all
positives, a capped dark-negative query universe, and a phylum-stratified
characterised-at-T0 reference. The raw-input budget any novelty scorer
needs, independent of whether Track 1 picks EVT or density-based scoring.

Usage:
    uv run python scripts/build_embedding_sample.py --release R207
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from darkmatter.embedding_sample import build_embedding_sample

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207")
    parser.add_argument("--dark-negative-cap", type=int, default=30000)
    parser.add_argument("--characterised-cap", type=int, default=30000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    labels = pd.read_csv(PROC_ROOT / f"gtdb_{args.release}" / "panel_protein_labels.csv")
    panel = pd.read_csv(PROC_ROOT / f"gtdb_{args.release}" / "genome_panel.csv")

    sample = build_embedding_sample(
        labels,
        panel,
        dark_negative_cap=args.dark_negative_cap,
        characterised_cap=args.characterised_cap,
        seed=args.seed,
    )

    out_path = PROC_ROOT / f"gtdb_{args.release}" / "embedding_sample.csv"
    sample.to_csv(out_path, index=False)

    print(f"{len(sample)} sequences sampled")
    print(sample["category"].value_counts().to_string())
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
