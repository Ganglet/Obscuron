#!/usr/bin/env python
"""CLI: Precision@K evaluation for the Layer 1 novelty scorer (P2-D3),
against Track 1's frozen config/scorer.yaml. Reports lift over the
scored set's own positive base rate for both the calibrated headline
and the raw-distance baseline (P2-D1) it must beat.

Usage:
    uv run python scripts/build_protein_table.py --release R207 --model esm2
    uv run python scripts/evaluate_scorer.py --arm esm2 --config config/scorer.yaml
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from darkmatter.experiment_log import log_experiment
from darkmatter.scoring.gpd import diagnostics
from darkmatter.scoring.knn import l2_normalize
from darkmatter.scoring.scorer import NoveltyScorer

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"
RESULTS_ROOT = Path(__file__).resolve().parents[1] / "results"


def precision_at_k(scores: np.ndarray, labels: np.ndarray, k: int) -> float:
    """Larger score = more novel = ranked first."""
    order = np.argsort(-np.asarray(scores), kind="stable")
    k = min(k, len(order))
    return float(np.asarray(labels)[order[:k]].mean())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", choices=["esm2", "genos_m"], default="esm2")
    parser.add_argument("--config", type=Path, default=Path("config/scorer.yaml"))
    parser.add_argument("--release", default="R207")
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    arm_cfg = config["arms"][args.arm]
    arm_slug = args.arm.replace("_", "-")

    embeddings = np.load(arm_cfg["embeddings"]).astype(np.float64)
    embeddings = l2_normalize(embeddings)

    table_path = PROC_ROOT / f"gtdb_{args.release}" / f"{arm_slug}_protein_table.csv"
    if not table_path.exists():
        raise SystemExit(f"{table_path} not found -- run scripts/build_protein_table.py first")
    table = pd.read_csv(table_path)
    assert len(table) == len(embeddings), "protein table / embedding matrix row count mismatch"

    ref_mask = (table["role"] == "reference").to_numpy()
    query_mask = (table["role"] == "query").to_numpy()
    reference = embeddings[ref_mask]
    query = embeddings[query_mask]
    query_table = table[query_mask].reset_index(drop=True)

    print(f"{len(reference)} reference, {len(query)} query", flush=True)

    scorer = NoveltyScorer(config)
    print(f"fitting scorer (k={scorer.k}, threshold_quantile={scorer.threshold_quantile})...", flush=True)
    scorer.fit(reference)
    tail = scorer.tail
    print(
        f"  u={tail.u:.4f} xi={tail.xi:.4f} beta={tail.beta:.4f} "
        f"n_exceedances={tail.n_exceedances}/{tail.n_calibration}",
        flush=True,
    )

    print("scoring queries (headline + baseline)...", flush=True)
    scored = scorer.score(query)
    baseline_distance = scorer.score_baseline(query)

    is_positive = query_table["is_positive"].to_numpy()
    set_positive_rate = float(is_positive.mean())
    print(
        f"set positive rate: {set_positive_rate:.4f} "
        "(natural dark-population density ~1.4%; this set is enriched by the P1-D6 embedding-budget cap)",
        flush=True,
    )

    eval_cfg = config["evaluation"]["precision_at_k"]
    precision_curve = {k: precision_at_k(scored["novelty"].to_numpy(), is_positive, k) for k in eval_cfg["K"]}
    baseline_curve = {k: precision_at_k(baseline_distance, is_positive, k) for k in eval_cfg["K"]}
    lift_curve = {k: (p / set_positive_rate if set_positive_rate > 0 else float("nan")) for k, p in precision_curve.items()}
    baseline_lift_curve = {
        k: (p / set_positive_rate if set_positive_rate > 0 else float("nan")) for k, p in baseline_curve.items()
    }

    print("Precision@K -- headline (GPD-calibrated ranking):", flush=True)
    for k in eval_cfg["K"]:
        print(f"  P@{k} = {precision_curve[k]:.4f} (lift {lift_curve[k]:.2f}x)", flush=True)
    print("Precision@K -- baseline (raw kNN distance ranking):", flush=True)
    for k in eval_cfg["K"]:
        print(f"  P@{k} = {baseline_curve[k]:.4f} (lift {baseline_lift_curve[k]:.2f}x)", flush=True)

    diag = diagnostics(scorer.loo_distances, tail)

    RESULTS_ROOT.mkdir(exist_ok=True)
    ranked = query_table.copy()
    ranked["distance"] = scored["distance"].to_numpy()
    ranked["pvalue"] = scored["pvalue"].to_numpy()
    ranked["novelty"] = scored["novelty"].to_numpy()
    ranked["baseline_distance"] = baseline_distance
    ranked = ranked.sort_values("novelty", ascending=False)

    ranked_path = RESULTS_ROOT / f"scorer_{args.arm}_precisionK.csv"
    ranked.to_csv(ranked_path, index=False)

    summary = {
        "arm": args.arm,
        "n_reference": len(reference),
        "n_query": len(query),
        "n_positive": int(is_positive.sum()),
        "set_positive_rate": set_positive_rate,
        "gpd": {
            "u": tail.u, "xi": tail.xi, "beta": tail.beta,
            "zeta_u": tail.zeta_u, "n_exceedances": tail.n_exceedances, "n_calibration": tail.n_calibration,
        },
        "gpd_goodness_of_fit_ks_pvalue": diag["ks_pvalue"],
        "precision_at_k": precision_curve,
        "lift_at_k": lift_curve,
        "baseline_precision_at_k": baseline_curve,
        "baseline_lift_at_k": baseline_lift_curve,
    }
    summary_path = RESULTS_ROOT / f"scorer_{args.arm}_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"wrote {ranked_path}", flush=True)
    print(f"wrote {summary_path}", flush=True)

    log_experiment(
        title=f"scorer evaluation -- precision@k ({args.arm})",
        config=f"k={scorer.k}, threshold_quantile={scorer.threshold_quantile}, K={eval_cfg['K']} (frozen, config/scorer.yaml)",
        result=(
            f"P@K={precision_curve}, lift={lift_curve}; baseline P@K={baseline_curve}, "
            f"baseline lift={baseline_lift_curve}; set_positive_rate={set_positive_rate:.4f}; "
            f"gpd ks_pvalue={diag['ks_pvalue']:.4f}"
        ),
        next_step="held-out-family AUROC + calibration via scripts/heldout_family_eval.py",
    )


if __name__ == "__main__":
    main()
