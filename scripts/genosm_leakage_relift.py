#!/usr/bin/env python
"""Re-test the Layer-4 Genos-m convergence axis's positive-rate enrichment
under the leakage-tightened 36->37 boundary (P4-D11), instead of the main
35->37 boundary used everywhere else in the project.

Background: the Genos-m axis's top-decile novelty enriches for positives
(lift 1.30x, P4-D8) -- the only axis in the project that points this
direction -- plausibly because Genos-m's pretraining (GTDB R220, Apr 2024)
sits inside the main 35.0->37.0 (Nov 2021 -> Jun 2024) boundary window.
This script asks: does the enrichment persist, shrink, or vanish when
"positive" is redefined under the tightest boundary this project's
Pfam releases can support (36.0/Sep 2023 -> 37.0/Jun 2024), now computed
at full panel scale (panel_protein_labels_36_37.py)?

Usage:
    uv run python scripts/genosm_leakage_relift.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "gtdb_R207"

GENOSM_SCORES = ROOT / "results" / "genos-m_novelty_scores.csv"
LABELS_36_37 = PROC / "panel_protein_labels_36_37.csv"
LABELS_35_37 = PROC / "panel_protein_labels.csv"
OUT_JSON = PROC / "genosm_leakage_relift.json"

QUANTILE = 0.90  # matches config/convergence.yaml's top-decile convention


def lift(df: pd.DataFrame, positive_col: str, score_col: str) -> dict:
    base_rate = df[positive_col].mean()
    thresh = df[score_col].quantile(QUANTILE)
    top = df[df[score_col] >= thresh]
    top_rate = top[positive_col].mean()
    return {
        "n_total": int(len(df)),
        "n_positive": int(df[positive_col].sum()),
        "base_rate": float(base_rate),
        "n_top_decile": int(len(top)),
        "top_decile_positive_rate": float(top_rate),
        "lift": float(top_rate / base_rate) if base_rate > 0 else float("nan"),
    }


def main() -> None:
    scores = pd.read_csv(GENOSM_SCORES)
    print(f"{len(scores)} genos-m-scored queries", flush=True)

    labels_36_37 = pd.read_csv(LABELS_36_37, usecols=["protein_id", "positive_36_37"])
    labels_35_37 = pd.read_csv(LABELS_35_37, usecols=["protein_id", "positive_proxy"])

    merged = scores.merge(labels_36_37, on="protein_id", how="left")
    merged = merged.merge(labels_35_37, on="protein_id", how="left")
    n_missing_36_37 = merged["positive_36_37"].isna().sum()
    n_missing_35_37 = merged["positive_proxy"].isna().sum()
    print(f"{n_missing_36_37} queries missing a 36-37 label, {n_missing_35_37} missing a 35-37 label (dropped)", flush=True)
    merged = merged.dropna(subset=["positive_36_37", "positive_proxy"])
    merged["positive_36_37"] = merged["positive_36_37"].astype(bool)
    merged["positive_proxy"] = merged["positive_proxy"].astype(bool)

    result_main = lift(merged, "positive_proxy", "genos_m_novelty")
    result_tight = lift(merged, "positive_36_37", "genos_m_novelty")

    print(f"\n=== main boundary (35->37), same query set, recomputed for consistency ===")
    print(json.dumps(result_main, indent=2))
    print(f"\n=== leakage-tightened boundary (36->37) ===")
    print(json.dumps(result_tight, indent=2))

    n_overlap = int((merged["positive_36_37"] & merged["positive_proxy"]).sum())
    print(f"\n{n_overlap} queries are positive under BOTH boundaries (of {int(merged['positive_36_37'].sum())} positive-36-37)")

    out = {
        "n_queries_matched": int(len(merged)),
        "main_boundary_35_37": result_main,
        "leakage_tightened_36_37": result_tight,
        "n_positive_both_boundaries": n_overlap,
    }
    OUT_JSON.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {OUT_JSON}")


if __name__ == "__main__":
    main()
