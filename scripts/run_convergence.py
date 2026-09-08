#!/usr/bin/env python
"""Layer 4 driver (Phase-4 extension): multi-signal convergence over independent
novelty axes (SS4 principle: confident novelty only where independent lines agree).

Axis 1 (always): Layer-1 ESM-2 EVT embedding novelty (esm2_novelty_scores.csv).
Axis 2 (always): genomic-context novelty (contig gene-order; embedding-free).
Axis 3+ (optional, --extra-axis): any per-query novelty CSV (protein_id + a score
    column) -- e.g. composition_novelty.csv now, or a cloud-produced
    genos-m_novelty_scores.csv later. Same wiring for both.

Reports, on the shared query set:
  1. Pairwise independence of the axes (Spearman + top-K Jaccard) -- the
     precondition for convergence to be real multi-evidence.
  2. The high-confidence convergent set (top-quantile on ALL axes): its size and
     positive-rate vs each single axis and vs the base rate.

Frozen params in config/convergence.yaml. Usage:
    uv run python scripts/run_convergence.py
    uv run python scripts/run_convergence.py --extra-axis data/processed/gtdb_R207/composition_novelty.csv
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from darkmatter.convergence import (
    convergence_stats,
    genomic_context_novelty,
    high_confidence_convergent_multi,
)

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "gtdb_R207"
SCORES = PROC / "esm2_novelty_scores.csv"
LABELS = PROC / "panel_protein_labels.csv"
OUT = PROC / "convergence_results.json"


def _rate(mask: np.ndarray, is_pos: np.ndarray) -> tuple[int, float]:
    n = int(mask.sum())
    return n, (float(is_pos[mask].mean()) if n else float("nan"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("--min-neighbors", type=int, default=4)
    ap.add_argument("--quantile", type=float, default=0.90)
    ap.add_argument("--top-ks", default="50,100,500,1000")
    ap.add_argument("--extra-axis", default=None,
                    help="CSV (protein_id + one score column) to add as a 3rd+ axis")
    ap.add_argument("--extra-name", default=None, help="label for the extra axis")
    args = ap.parse_args()
    top_ks = [int(x) for x in args.top_ks.split(",")]

    scores = pd.read_csv(SCORES)
    print(f"{len(scores)} scored dark queries (EVT novelty)", flush=True)
    labels = pd.read_csv(LABELS, usecols=["protein_id", "dark_at_t0"])
    context = genomic_context_novelty(labels, args.window, args.min_neighbors)
    print(f"context novelty defined for {len(context)} genes", flush=True)
    scores["context_novelty"] = scores["protein_id"].map(context)

    axis_cols = {"evt": "evt_novelty_score", "context": "context_novelty"}
    if args.extra_axis:
        extra = pd.read_csv(args.extra_axis)
        score_col = [c for c in extra.columns if c != "protein_id"][0]
        name = args.extra_name or score_col.replace("_novelty", "")
        scores = scores.merge(extra[["protein_id", score_col]], on="protein_id", how="left")
        axis_cols[name] = score_col
        print(f"extra axis '{name}' from {Path(args.extra_axis).name} ({score_col})", flush=True)

    merged = scores.dropna(subset=list(axis_cols.values())).reset_index(drop=True)
    is_pos = merged["is_positive"].astype(str).eq("True").to_numpy()
    axes = {name: merged[col].to_numpy() for name, col in axis_cols.items()}
    print(f"{len(merged)} queries with ALL {len(axes)} signals; {int(is_pos.sum())} positives", flush=True)

    # pairwise independence
    pairwise = {}
    for a, b in combinations(axes, 2):
        pairwise[f"{a}~{b}"] = convergence_stats(axes[a], axes[b], top_ks)

    # convergent set = high on ALL axes
    both = high_confidence_convergent_multi(list(axes.values()), args.quantile)
    base_n, base_rate = _rate(np.ones(len(merged), bool), is_pos)
    singles = {}
    for name, v in axes.items():
        n, r = _rate(v >= np.quantile(v, args.quantile), is_pos)
        singles[name] = {"n": n, "positive_rate": r, "lift": r / base_rate}
    conv_n, conv_rate = _rate(both, is_pos)

    result = {
        "params": vars(args),
        "axes": list(axes),
        "n_queries_all_signals": int(len(merged)),
        "base_positive_rate": base_rate,
        "pairwise_independence": pairwise,
        "single_axis_top": singles,
        "convergent_all": {"n": conv_n, "positive_rate": conv_rate,
                            "lift": (conv_rate / base_rate) if conv_rate == conv_rate else None},
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print(f"\n=== pairwise axis independence (near 0 Spearman = independent = real multi-evidence) ===")
    for pair, s in pairwise.items():
        js = "  ".join(f"K{k}={s['top_jaccard'][k]:.3f}" for k in top_ks)
        print(f"{pair:24} Spearman {s['spearman']:+.4f}   Jaccard {js}")

    print(f"\n=== what agreement buys (positive-rate; base = near-known bias) ===")
    print(f"base rate               : {base_rate:.4f} ({base_n})")
    for name, s in singles.items():
        print(f"{name+'-top '+str(int(100*args.quantile))+'%':<24}: {s['positive_rate']:.4f}  lift {s['lift']:.2f}x ({s['n']})")
    print(f"{'CONVERGENT (all '+str(len(axes))+')':<24}: {conv_rate:.4f}  lift {conv_rate/base_rate:.2f}x ({conv_n})")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
