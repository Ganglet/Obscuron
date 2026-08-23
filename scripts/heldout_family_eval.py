#!/usr/bin/env python
"""CLI: held-out-family AUROC + calibration for the Layer 1 novelty
scorer (P2-D4), against Track 1's frozen config/scorer.yaml.

Withholds one Pfam-35 family at a time from the reference, refits the
scorer on the reduced reference, and scores the withheld family's
members (label "novel") against a matched sample of retained knowns
(label "known") -- one control per novel member, nearest length_aa
within the same phylum, sampled without replacement (P2-D4's
length/phylum matching; family-size matching is implicit in drawing
exactly as many controls as novel members per family). Aggregates into
an AUROC distribution and a pooled reliability diagram -- the controlled
ground truth calibration is certified on, never the selection-biased
positives (P2-D3).

Usage:
    uv run python scripts/build_protein_table.py --release R207 --model esm2
    uv run python scripts/heldout_family_eval.py --arm esm2 --config config/scorer.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats

from darkmatter.experiment_log import log_experiment
from darkmatter.scoring.knn import l2_normalize
from darkmatter.scoring.scorer import NoveltyScorer

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"
RESULTS_ROOT = Path(__file__).resolve().parents[1] / "results"


def _auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=bool)
    n_pos, n_neg = int(labels.sum()), int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = stats.rankdata(scores)
    sum_ranks_pos = ranks[labels].sum()
    return float((sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def _match_controls(family_rows: pd.DataFrame, pool: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """One control per novel member: nearest length_aa within the same
    phylum, sampled without replacement from the retained pool. rng is
    accepted for interface consistency (deterministic nearest-match here,
    no randomness needed once phylum-filtered)."""
    del rng
    chosen_idx: list[int] = []
    used: set[int] = set()
    for _, row in family_rows.iterrows():
        candidates = pool[(pool["phylum"] == row["phylum"]) & (~pool.index.isin(used))]
        if candidates.empty:
            candidates = pool[~pool.index.isin(used)]
        if candidates.empty:
            continue
        diffs = (candidates["length_aa"] - row["length_aa"]).abs()
        best_idx = diffs.idxmin()
        chosen_idx.append(best_idx)
        used.add(best_idx)
    return pool.loc[chosen_idx]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", choices=["esm2", "genos_m"], default="esm2")
    parser.add_argument("--config", type=Path, default=Path("config/scorer.yaml"))
    parser.add_argument("--release", default="R207")
    parser.add_argument(
        "--max-families", type=int, default=None, help="override config's max_families, for a quick test run"
    )
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    arm_cfg = config["arms"][args.arm]
    arm_slug = args.arm.replace("_", "-")
    ho_cfg = config["evaluation"]["heldout_family"]
    seed = config["seed"]

    embeddings = np.load(arm_cfg["embeddings"]).astype(np.float64)
    embeddings = l2_normalize(embeddings)

    table_path = PROC_ROOT / f"gtdb_{args.release}" / f"{arm_slug}_protein_table.csv"
    if not table_path.exists():
        raise SystemExit(f"{table_path} not found -- run scripts/build_protein_table.py first")
    table = pd.read_csv(table_path)
    assert len(table) == len(embeddings), "protein table / embedding matrix row count mismatch"

    ref_mask = (table["role"] == "reference").to_numpy()
    ref_table = table[ref_mask].reset_index(drop=True)
    ref_table["_emb_idx"] = np.where(ref_mask)[0]

    family_sizes = ref_table.dropna(subset=["pfam_family"]).groupby("pfam_family").size()
    candidate_families = family_sizes[family_sizes >= ho_cfg["min_members"]].index.tolist()
    print(f"{len(candidate_families)} candidate families with >= {ho_cfg['min_members']} members", flush=True)

    max_families = args.max_families if args.max_families is not None else ho_cfg["max_families"]
    rng = np.random.default_rng(seed)
    if len(candidate_families) > max_families:
        strat_df = pd.DataFrame({"pfam_family": candidate_families})
        strat_df["family_size"] = strat_df["pfam_family"].map(family_sizes)
        n_buckets = min(5, strat_df["family_size"].nunique())
        strat_df["size_bucket"] = pd.qcut(strat_df["family_size"], q=n_buckets, duplicates="drop")
        majority_phylum = ref_table.dropna(subset=["pfam_family"]).groupby("pfam_family")["phylum"].agg(
            lambda s: s.mode().iloc[0]
        )
        strat_df["phylum"] = strat_df["pfam_family"].map(majority_phylum)
        frac = max_families / len(strat_df)
        sampled = strat_df.groupby(["phylum", "size_bucket"], group_keys=False, observed=True).apply(
            lambda g: g.sample(n=max(1, round(len(g) * frac)), random_state=rng)
        )
        selected_families = sampled["pfam_family"].tolist()[:max_families]
    else:
        selected_families = candidate_families
    print(f"evaluating {len(selected_families)} held-out families", flush=True)

    per_family_results = []
    pooled_pvalues: list[float] = []
    pooled_labels: list[int] = []

    for i, family in enumerate(selected_families, 1):
        family_rows = ref_table[ref_table["pfam_family"] == family]
        retained = ref_table[ref_table["pfam_family"] != family]

        reduced_reference = embeddings[retained["_emb_idx"].to_numpy()]
        novel_embeddings = embeddings[family_rows["_emb_idx"].to_numpy()]

        controls = _match_controls(family_rows, retained, rng)
        known_embeddings = embeddings[controls["_emb_idx"].to_numpy()]

        scorer = NoveltyScorer(config)
        scorer.fit(reduced_reference)
        novel_scored = scorer.score(novel_embeddings)
        known_scored = scorer.score(known_embeddings)

        scores = np.concatenate([novel_scored["novelty"].to_numpy(), known_scored["novelty"].to_numpy()])
        labels = np.concatenate([np.ones(len(novel_scored)), np.zeros(len(known_scored))])
        family_auroc = _auroc(scores, labels)

        per_family_results.append(
            {"pfam_family": family, "n_novel": len(family_rows), "n_known": len(controls), "auroc": family_auroc}
        )
        pooled_pvalues.extend(novel_scored["pvalue"].tolist())
        pooled_pvalues.extend(known_scored["pvalue"].tolist())
        pooled_labels.extend([1] * len(novel_scored))
        pooled_labels.extend([0] * len(known_scored))

        if i % 25 == 0 or i == len(selected_families):
            print(f"  [{i}/{len(selected_families)}] families done", flush=True)

    results_df = pd.DataFrame(per_family_results)
    RESULTS_ROOT.mkdir(exist_ok=True)
    auroc_path = RESULTS_ROOT / f"heldout_{args.arm}_auroc.csv"
    results_df.to_csv(auroc_path, index=False)

    median_auroc = float(results_df["auroc"].median())
    mean_auroc = float(results_df["auroc"].mean())
    print(f"held-out-family AUROC: median={median_auroc:.4f}, mean={mean_auroc:.4f}, n_families={len(results_df)}", flush=True)

    pooled = pd.DataFrame({"pvalue": pooled_pvalues, "label": pooled_labels})
    pooled["novelty"] = 1.0 - pooled["pvalue"]
    n_bins = config["evaluation"]["calibration"]["bins"]
    bins = np.linspace(0, 1, n_bins + 1)
    pooled["bin"] = pd.cut(pooled["novelty"], bins=bins, include_lowest=True)
    calibration = pooled.groupby("bin", observed=True)["label"].agg(["mean", "count"])
    print("calibration (held-out-family ground truth):", flush=True)
    print(calibration.to_string(), flush=True)

    calibration_path = RESULTS_ROOT / f"heldout_{args.arm}_calibration.csv"
    calibration.to_csv(calibration_path)

    fig_dir = RESULTS_ROOT / "figures"
    fig_dir.mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.hist(results_df["auroc"].dropna(), bins=30)
    ax.axvline(0.5, color="k", linestyle="--", alpha=0.5, label="chance")
    ax.set_xlabel("AUROC")
    ax.set_ylabel("count of held-out families")
    ax.set_title(f"Held-out-family AUROC distribution ({args.arm}, n={len(results_df)})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.arm}_heldout_auroc_distribution.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    bin_centers = calibration.index.map(lambda iv: iv.mid).astype(float)
    ax.plot(bin_centers, calibration["mean"], "o-", label="observed novel-fraction")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.4, label="perfect calibration")
    ax.set_xlabel("predicted novelty score (bin center)")
    ax.set_ylabel("observed novel-fraction (held-out-family ground truth)")
    ax.set_title(f"Calibration -- held-out-family ground truth ({args.arm})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.arm}_heldout_calibration.png", dpi=150)
    plt.close(fig)

    print(f"wrote {auroc_path}", flush=True)
    print(f"wrote {calibration_path}", flush=True)
    print(f"wrote figures to {fig_dir}", flush=True)

    log_experiment(
        title=f"held-out-family evaluation ({args.arm}) -- AUROC distribution + calibration",
        config=f"min_members={ho_cfg['min_members']}, max_families={max_families}, seed={seed} (frozen, config/scorer.yaml)",
        result=f"{len(results_df)} families: median AUROC={median_auroc:.4f}, mean={mean_auroc:.4f}; calibration in {calibration_path}",
        next_step="report to Track 1 for interpretation against the acceptance sanity checks in docs/Track2_Phase2_scoring_handoff.md",
    )


if __name__ == "__main__":
    main()
