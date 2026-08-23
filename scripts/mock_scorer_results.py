#!/usr/bin/env python
"""Generate MOCK schema-correct scorer results to smoke-test the analysis harness.

NOT a scorer, and never used for any claim. Fabricates a plausible *healthy* result
matching the frozen result schema (docs/Track2_Phase2_scoring_handoff.md) so
scripts/analyze_scorer_results.py can be built and verified before Track 2's real run
exists. Written to results/ (gitignored). Delete once real results land.

    uv run python scripts/mock_scorer_results.py --arm esm2
"""
from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PHYLA = [
    "Pseudomonadota", "Bacillota", "Actinomycetota",
    "Bacteroidota", "Patescibacteria", "Thermoproteota",
]


def main() -> None:
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("--arm", default="esm2")
    ap.add_argument("--out-dir", default=str(ROOT / "results"))
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    odir = Path(args.out_dir)
    odir.mkdir(parents=True, exist_ok=True)
    arm = args.arm

    # --- ranked queries: novelty skewed low; positives enriched at high novelty ---
    n = 25_000
    novelty = rng.beta(2.0, 5.0, n)
    is_pos = rng.random(n) < (0.02 + 0.30 * novelty)          # monotone signal -> lift > 1
    distance = np.clip(novelty + rng.normal(0, 0.35, n), 0, None)  # baseline: noisier proxy -> headline wins
    ranked = pd.DataFrame({
        "protein_id": [f"MOCK_{i:06d}" for i in range(n)],
        "novelty": novelty,
        "pvalue": 1.0 - novelty,
        "distance": distance,
        "is_positive": is_pos.astype(int),
        "phylum": rng.choice(PHYLA, n),
        "length_aa": rng.integers(80, 1024, n),
    }).sort_values("novelty", ascending=False).reset_index(drop=True)
    ranked.to_csv(odir / f"scorer_{arm}_ranked.csv", index=False)

    # --- GPD tail QQ: near-diagonal (good fit) ---
    q = np.linspace(0.02, 0.99, 40)
    theoretical = -np.log(1 - q)                              # unit-exponential-ish tail
    empirical = theoretical + rng.normal(0, 0.05, theoretical.size)
    pd.DataFrame({"empirical": empirical, "theoretical": theoretical}).to_csv(
        odir / f"scorer_{arm}_gpd_qq.csv", index=False)

    # --- held-out-family AUROC: median ~0.8 ---
    nf = 300
    auroc = np.clip(rng.normal(0.80, 0.09, nf), 0.5, 1.0)
    pd.DataFrame({
        "family_id": [f"PF{rng.integers(0, 99999):05d}" for _ in range(nf)],
        "n_members": rng.integers(5, 200, nf),
        "phylum": rng.choice(PHYLA, nf),
        "auroc": auroc,
    }).to_csv(odir / f"heldout_{arm}_auroc.csv", index=False)

    # --- calibration bins: near-diagonal ---
    pred = np.linspace(0.05, 0.95, 10)
    obs = np.clip(pred + rng.normal(0, 0.03, pred.size), 0, 1)
    pd.DataFrame({
        "predicted": pred,
        "observed": obs,
        "n": rng.integers(50, 400, pred.size),
    }).to_csv(odir / f"heldout_{arm}_calibration.csv", index=False)

    print(f"wrote MOCK results for arm={arm} to {odir}")


if __name__ == "__main__":
    main()
