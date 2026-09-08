#!/usr/bin/env python
"""Render the Phase-2 scorer figures + a summary table from Track 2's results.

Track 1 analysis: reads the frozen result schema (docs/Track2_Phase2_scoring_handoff.md),
writes the paper figures to figures/ and a markdown summary. Runs against real Track-2
output (the frozen result schema).

    uv run python scripts/analyze_scorer_results.py --arm esm2
"""
from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from darkmatter.analysis import figures as F

ROOT = Path(__file__).resolve().parents[1]
K_LIST = [50, 100, 500, 1000]


def main() -> None:
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("--arm", default="esm2")
    ap.add_argument("--results-dir", default=str(ROOT / "results"))
    ap.add_argument("--out-dir", default=str(ROOT / "figures"))
    args = ap.parse_args()

    rdir = Path(args.results_dir)
    odir = Path(args.out_dir)
    odir.mkdir(parents=True, exist_ok=True)
    arm = args.arm

    ranked = pd.read_csv(rdir / f"scorer_{arm}_ranked.csv")
    calib = pd.read_csv(rdir / f"heldout_{arm}_calibration.csv")
    auroc = pd.read_csv(rdir / f"heldout_{arm}_auroc.csv")
    qq = pd.read_csv(rdir / f"scorer_{arm}_gpd_qq.csv")

    panels = {
        f"scorer_{arm}_precision_at_k.png": F.precision_at_k_curve(ranked, K_LIST),
        f"scorer_{arm}_calibration.png": F.reliability_diagram(calib),
        f"scorer_{arm}_auroc.png": F.auroc_distribution(auroc),
        f"scorer_{arm}_gpd_qq.png": F.gpd_qq_plot(qq),
    }
    for name, fig in panels.items():
        fig.savefig(odir / name, bbox_inches="tight")
        print("wrote", odir / name)

    pk, base = F.precision_at_k(ranked, "novelty", K_LIST)
    bk, _ = F.precision_at_k(ranked, "distance", K_LIST)
    med = auroc["auroc"].median()
    lines = [
        f"# Phase 2 scorer summary ({arm})",
        "",
        f"- queries: {len(ranked)}   positives: {int(ranked['is_positive'].sum())}   "
        f"base rate: {base:.4f}",
        f"- held-out-family AUROC median: {med:.3f}   (families: {len(auroc)})",
        "",
        "| K | P@K novelty | lift | P@K distance | lift |",
        "|---|---|---|---|---|",
    ]
    for (_, h), (_, b) in zip(pk.iterrows(), bk.iterrows()):
        lines.append(
            f"| {int(h.K)} | {h.precision:.3f} | {h.lift:.2f} | {b.precision:.3f} | {b.lift:.2f} |"
        )
    (odir / f"scorer_{arm}_summary.md").write_text("\n".join(lines) + "\n")
    print("wrote", odir / f"scorer_{arm}_summary.md")


if __name__ == "__main__":
    main()
