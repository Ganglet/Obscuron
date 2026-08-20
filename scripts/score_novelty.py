#!/usr/bin/env python
"""CLI: Layer 1 calibrated novelty scorer (EVT on kNN cosine distance),
per Track 1's P2-D1/P2-D2 spec.

Reference = characterised_at_t0 (has a Pfam-35 GA hit); queries =
dark_negative + positive (no Pfam-35 hit at T0; positive = also
characterised by T1, the true label this scores against).

k (kNN neighbours) and the GPD threshold quantile are NOT yet locked by
Track 1 -- see docs/problems_and_decisions.md -- so these are
configurable flags with provisional defaults, not final evaluation
numbers. Precision@K, held-out-family AUROC, and the calibration check
still need Track 1 to fix K and the held-out family before real
evaluation can be reported (blueprint D5: metrics fixed before
implementation, never adjusted post hoc).

Usage:
    uv run python scripts/score_novelty.py --release R207 --model esm2
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from darkmatter.experiment_log import log_experiment
from darkmatter.novelty import (
    evt_novelty_score,
    fit_gpd_tail,
    knn_mean_distance,
    l2_normalize,
    leave_one_out_distances,
)

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207")
    parser.add_argument("--model", choices=["genos-m", "esm2"], default="esm2")
    parser.add_argument("--k", type=int, default=10, help="neighbors for kNN distance -- PROVISIONAL, not locked by Track 1 yet")
    parser.add_argument("--threshold-quantile", type=float, default=0.9, help="GPD tail threshold u -- PROVISIONAL")
    args = parser.parse_args()

    out_dir = PROC_ROOT / f"gtdb_{args.release}"
    embeddings = np.load(out_dir / f"{args.model}_panel_embeddings.npy").astype(np.float64)
    manifest = pd.read_csv(out_dir / f"{args.model}_panel_embeddings_manifest.csv")
    assert len(embeddings) == len(manifest), "embeddings/manifest row count mismatch"

    embeddings = l2_normalize(embeddings)

    ref_mask = (manifest["category"] == "characterised_at_t0").to_numpy()
    query_mask = manifest["category"].isin(["dark_negative", "positive"]).to_numpy()
    references = embeddings[ref_mask]
    queries = embeddings[query_mask]
    query_manifest = manifest[query_mask].reset_index(drop=True)

    print(f"{len(references)} reference (characterised-at-T0), {len(queries)} queries (dark-at-T0)", flush=True)

    print(f"leave-one-out calibration on reference (k={args.k})...", flush=True)
    loo_distances = leave_one_out_distances(references, k=args.k)

    print(f"fitting GPD tail (threshold quantile={args.threshold_quantile})...", flush=True)
    gpd = fit_gpd_tail(loo_distances, threshold_quantile=args.threshold_quantile)
    print(
        f"  u={gpd['u']:.4f}, shape={gpd['shape']:.4f}, scale={gpd['scale']:.4f}, "
        f"n_exceedances={gpd['n_exceedances']}/{gpd['n_calibration']}",
        flush=True,
    )

    print("scoring queries (raw kNN distance + EVT-calibrated)...", flush=True)
    query_distances = knn_mean_distance(queries, references, k=args.k)
    evt_scores = evt_novelty_score(query_distances, gpd)

    result = query_manifest.copy()
    result["raw_knn_distance"] = query_distances
    result["evt_novelty_score"] = evt_scores
    result["is_positive"] = result["category"] == "positive"

    out_path = out_dir / f"{args.model}_novelty_scores.csv"
    result.to_csv(out_path, index=False)
    print(f"wrote {out_path}", flush=True)

    n_positive = int(result["is_positive"].sum())
    positive_mean_score = float(result.loc[result["is_positive"], "evt_novelty_score"].mean())
    negative_mean_score = float(result.loc[~result["is_positive"], "evt_novelty_score"].mean())
    print(
        f"{len(result)} queries scored, {n_positive} true positives "
        f"(mean score {positive_mean_score:.4f}) vs {len(result) - n_positive} still-dark "
        f"(mean score {negative_mean_score:.4f})",
        flush=True,
    )

    log_experiment(
        title=f"layer 1 evt novelty scorer ({args.model}, {args.release} panel)",
        config=f"k={args.k}, gpd_threshold_quantile={args.threshold_quantile} (PROVISIONAL, not locked by Track 1)",
        result=(
            f"{len(references)} reference, {len(queries)} queries scored ({n_positive} true positives); "
            f"gpd u={gpd['u']:.4f} shape={gpd['shape']:.4f}; "
            f"mean novelty positive={positive_mean_score:.4f} vs still-dark={negative_mean_score:.4f}"
        ),
        next_step="Precision@K, held-out-family AUROC, and calibration check still need Track 1 to fix K and the held-out family before real evaluation numbers can be reported",
    )


if __name__ == "__main__":
    main()
