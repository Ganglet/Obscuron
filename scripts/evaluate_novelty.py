#!/usr/bin/env python
"""CLI: evaluate the Layer 1 EVT novelty scorer against blueprint D5's
three metrics -- Precision@K, held-out-family AUROC, and a calibration
check -- plus the score-distribution/ROC/calibration figures Track 2
owes for Phase 2.

K (for Precision@K), the held-out family, and the scorer's own k/u
(score_novelty.py) are NOT yet locked by Track 1 -- see
docs/problems_and_decisions.md. This runs with configurable, PROVISIONAL
defaults so the mechanism is ready to go the instant real values land;
treat none of these numbers as final per blueprint D5.

Usage:
    uv run python scripts/score_novelty.py --release R207 --model esm2
    uv run python scripts/evaluate_novelty.py --release R207 --model esm2
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyhmmer.easel import Alphabet

from darkmatter.data.hmmscan import load_protein_sequences_multi, scan_against_pfam
from darkmatter.experiment_log import log_experiment
from darkmatter.novelty import (
    auroc,
    evt_novelty_score,
    fit_gpd_tail,
    knn_mean_distance,
    l2_normalize,
    leave_one_out_distances,
    precision_at_k,
)

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"
RAW_ROOT = Path(__file__).resolve().parents[1] / "data" / "raw"
PFAM35_HMM = RAW_ROOT / "pfam_35.0" / "Pfam-A.hmm.gz"


def _pick_held_out_family(out_dir: Path, ref_protein_ids: set[str], min_family_size: int, cpus: int) -> tuple[str, set[str]]:
    """Scan the reference proteins against Pfam-35 (they all have >=1 hit
    by definition of characterised-at-T0) and pick the largest family of
    unambiguous single-hit members as the held-out family."""
    manifest = pd.read_csv(out_dir / "genome_panel.csv")
    proteins_dir = out_dir / "panel_proteins"

    # Only need to load genomes that actually contain a reference protein.
    ref_genomes = set()
    embed_manifest = pd.read_csv(out_dir / "esm2_panel_embeddings_manifest.csv")
    ref_rows = embed_manifest[embed_manifest["protein_id"].isin(ref_protein_ids)]
    ref_genomes = set(ref_rows["genome_accession"])

    faa_files = [proteins_dir / f"{acc}_protein.faa" for acc in ref_genomes]
    faa_files = [f for f in faa_files if f.exists()]

    alphabet = Alphabet.amino()
    sequences = load_protein_sequences_multi(faa_files, alphabet)
    # Restrict to just the reference proteins we actually need scored.
    ref_sequences = [s for s in sequences if s.name in ref_protein_ids]

    print(f"scanning {len(ref_sequences)} reference proteins against Pfam-35 for family labels...", flush=True)
    from pyhmmer.easel import DigitalSequenceBlock

    ref_block = DigitalSequenceBlock(alphabet, ref_sequences)
    hits = scan_against_pfam(ref_block, PFAM35_HMM, cpus=cpus)

    from collections import defaultdict

    by_family: dict[str, list[str]] = defaultdict(list)
    for protein_id, families in hits.items():
        if len(families) == 1:  # unambiguous
            by_family[families[0]].append(protein_id)

    candidates = {fam: ids for fam, ids in by_family.items() if len(ids) >= min_family_size}
    if not candidates:
        raise ValueError(f"no Pfam-35 family with >= {min_family_size} unambiguous reference members found")
    held_out_family = max(candidates, key=lambda f: len(candidates[f]))
    return held_out_family, set(candidates[held_out_family])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207")
    parser.add_argument("--model", choices=["genos-m", "esm2"], default="esm2")
    parser.add_argument("--k", type=int, default=10, help="kNN neighbors -- PROVISIONAL, must match score_novelty.py's --k")
    parser.add_argument("--threshold-quantile", type=float, default=0.9, help="GPD tail threshold u -- PROVISIONAL")
    parser.add_argument("--min-family-size", type=int, default=20, help="minimum unambiguous members to be eligible as the held-out family")
    parser.add_argument("--precision-k", type=int, nargs="+", default=[50, 100, 500, 1000, 5000], help="K values for Precision@K -- PROVISIONAL, none locked by Track 1 yet")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--calibration-bins", type=int, default=10)
    args = parser.parse_args()

    out_dir = PROC_ROOT / f"gtdb_{args.release}"

    # ---- Precision@K (reuses the already-scored queries) ----
    scores_path = out_dir / f"{args.model}_novelty_scores.csv"
    if not scores_path.exists():
        raise SystemExit(f"{scores_path} not found -- run scripts/score_novelty.py first")
    scored = pd.read_csv(scores_path)
    print(f"loaded {len(scored)} scored queries from {scores_path}", flush=True)

    precision_results = {
        k: precision_at_k(scored["evt_novelty_score"].to_numpy(), scored["is_positive"].to_numpy(), k)
        for k in args.precision_k
    }
    baseline_rate = float(scored["is_positive"].mean())
    print(f"Precision@K (baseline positive rate = {baseline_rate:.4f}):", flush=True)
    for k, p in precision_results.items():
        print(f"  P@{k} = {p:.4f} ({p / baseline_rate:.2f}x baseline)", flush=True)

    # ---- Held-out-family AUROC ----
    embeddings = np.load(out_dir / f"{args.model}_panel_embeddings.npy").astype(np.float64)
    manifest = pd.read_csv(out_dir / f"{args.model}_panel_embeddings_manifest.csv")
    embeddings = l2_normalize(embeddings)

    ref_mask = (manifest["category"] == "characterised_at_t0").to_numpy()
    ref_protein_ids = set(manifest.loc[ref_mask, "protein_id"])

    held_out_family, held_out_ids = _pick_held_out_family(out_dir, ref_protein_ids, args.min_family_size, args.cpus)
    print(f"held-out family: {held_out_family} ({len(held_out_ids)} unambiguous members)", flush=True)

    is_held_out = manifest["protein_id"].isin(held_out_ids).to_numpy()
    reduced_ref_mask = ref_mask & ~is_held_out
    references = embeddings[reduced_ref_mask]
    held_out_embeddings = embeddings[ref_mask & is_held_out]

    print(f"recalibrating on reduced reference ({len(references)}, family removed)...", flush=True)
    loo = leave_one_out_distances(references, k=args.k)
    gpd = fit_gpd_tail(loo, threshold_quantile=args.threshold_quantile)

    # "Normal" comparison class: an equal-sized random sample of the
    # reduced reference's own members (should NOT be flagged anomalous).
    rng = np.random.default_rng(42)
    normal_sample_idx = rng.choice(len(references), size=min(len(held_out_embeddings), len(references)), replace=False)
    normal_sample = references[normal_sample_idx]

    held_out_distances = knn_mean_distance(held_out_embeddings, references, k=args.k)
    normal_distances = knn_mean_distance(normal_sample, references, k=args.k)
    held_out_scores = evt_novelty_score(held_out_distances, gpd)
    normal_scores = evt_novelty_score(normal_distances, gpd)

    auroc_scores = np.concatenate([held_out_scores, normal_scores])
    auroc_labels = np.concatenate([np.ones(len(held_out_scores)), np.zeros(len(normal_scores))])
    held_out_auroc = auroc(auroc_scores, auroc_labels)
    print(f"held-out-family AUROC: {held_out_auroc:.4f} (0.5 = chance, 1.0 = perfect separation)", flush=True)

    # ---- Calibration check ----
    bins = np.linspace(0, 1, args.calibration_bins + 1)
    scored["score_bin"] = pd.cut(scored["evt_novelty_score"], bins=bins, include_lowest=True)
    calibration = scored.groupby("score_bin", observed=True)["is_positive"].agg(["mean", "count"])
    print("calibration (predicted-score bin vs observed positive rate):", flush=True)
    print(calibration.to_string(), flush=True)

    # ---- Diagnostics: figures ----
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(scored.loc[~scored["is_positive"], "evt_novelty_score"], bins=30, alpha=0.6, label="still-dark", density=True)
    ax.hist(scored.loc[scored["is_positive"], "evt_novelty_score"], bins=30, alpha=0.6, label="positive (characterised-by-T1)", density=True)
    ax.set_xlabel("EVT novelty score")
    ax.set_ylabel("density")
    ax.set_title(f"Novelty score distribution ({args.model})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.model}_score_distribution.png", dpi=150)
    plt.close(fig)

    thresholds = np.linspace(0, 1, 200)
    tpr = [(held_out_scores >= t).mean() for t in thresholds]
    fpr = [(normal_scores >= t).mean() for t in thresholds]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(fpr, tpr, label=f"held-out family (AUROC={held_out_auroc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.4, label="chance")
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    ax.set_title(f"Held-out-family ROC ({held_out_family})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.model}_held_out_family_roc.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    bin_centers = calibration.index.map(lambda iv: iv.mid).astype(float)
    ax.plot(bin_centers, calibration["mean"], "o-", label="observed positive rate")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.4, label="perfect calibration")
    ax.set_xlabel("predicted novelty score (bin center)")
    ax.set_ylabel("observed positive rate")
    ax.set_title(f"Calibration check ({args.model})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / f"{args.model}_calibration.png", dpi=150)
    plt.close(fig)

    print(f"wrote figures to {fig_dir}", flush=True)

    log_experiment(
        title=f"layer 1 evaluation ({args.model}, {args.release} panel) -- Precision@K, held-out-family AUROC, calibration",
        config=(
            f"k={args.k}, gpd_threshold_quantile={args.threshold_quantile}, "
            f"held_out_family={held_out_family} ({len(held_out_ids)} members), "
            f"precision_k={args.precision_k} (ALL PROVISIONAL, not locked by Track 1)"
        ),
        result=(
            f"Precision@K: {precision_results} (baseline rate {baseline_rate:.4f}); "
            f"held-out-family AUROC={held_out_auroc:.4f}; "
            f"figures in {fig_dir}"
        ),
        next_step="waiting on Track 1 to lock k, GPD threshold, K, and the held-out family choice before these numbers are reportable as final",
    )


if __name__ == "__main__":
    main()
