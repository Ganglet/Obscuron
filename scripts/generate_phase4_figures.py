#!/usr/bin/env python
"""Figures for the three Phase-4 results that landed without one: Layer 4
(multi-signal convergence), Layer 5 (coding-structure statistics), and Layer 2
(ProstT5). The README already carries the numbers as tables (P4-D2/D4, P4-D5,
P4-D6) -- this just gives them the same visual treatment as the other four
results (matches src/darkmatter/analysis/figures.py's palette/style).

Layer 4 and Layer 5 plot real committed data (convergence_results.json,
statistical_results.json, coding_structure.csv). Layer 2 has no committed raw
file (ProstT5's .npz lives only on Angshuman's machine, gitignored per P1-D4) --
that panel reproduces the exact numbers already committed in the decision log /
README (P4-D6), not a re-derivation.

    uv run python scripts/generate_phase4_figures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "gtdb_R207"
OUT_DIRS = [PROC / "figures", ROOT / "results" / "figures"]
for d in OUT_DIRS:
    d.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
})
HEADLINE = "#3B6FB0"
BASELINE = "#C0603A"
THIRD = "#5B8C5A"
REF = "#8A8A8A"


def save(fig, name: str) -> None:
    for d in OUT_DIRS:
        fig.savefig(d / name)
    print(f"wrote {name} -> {[str(d) for d in OUT_DIRS]}")
    plt.close(fig)


def layer4_convergence() -> None:
    d = json.load(open(PROC / "convergence_results.json"))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.8), constrained_layout=True)

    pairs = list(d["pairwise_independence"].items())
    labels = [p[0].replace("~", " ~ ") for p in pairs]
    rhos = [p[1]["spearman"] for p in pairs]
    colors = [HEADLINE if abs(r) < 0.1 else BASELINE for r in rhos]
    ax1.barh(labels, rhos, color=colors)
    ax1.axvline(0, color=REF, lw=1)
    ax1.set_xlabel("Spearman rho")
    ax1.set_title("Pairwise axis independence")
    ax1.set_xlim(-0.1, 0.4)

    axes_order = ["evt", "context", "composition"]
    base = d["base_positive_rate"]
    lifts = [d["single_axis_top"][a]["lift"] for a in axes_order]
    ns = [d["single_axis_top"][a]["n"] for a in axes_order]
    conv_lift = d["convergent_all"]["lift"]
    conv_n = d["convergent_all"]["n"]
    labels2 = [f"{a}\n(n={n})" for a, n in zip(axes_order, ns)] + [f"3-way\n(n={conv_n})"]
    vals = lifts + [conv_lift]
    colors2 = [HEADLINE] * 3 + [THIRD]
    ax2.bar(labels2, vals, color=colors2)
    ax2.axhline(1.0, color=REF, ls=":", lw=1, label="base rate")
    ax2.set_ylabel(f"positive-rate lift (base {base:.3f})")
    ax2.set_title("Every axis depletes near-known positives")
    ax2.tick_params(axis="x", labelsize=8)
    ax2.legend(frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1.0, 0.93))
    fig.suptitle("Layer 4 — multi-signal convergence (P4-D2/D4)", fontsize=11)
    save(fig, "layer4_convergence_summary.png")


def layer5_coding_structure() -> None:
    d = json.load(open(PROC / "statistical_results.json"))
    cs = pd.read_csv(PROC / "coding_structure.csv")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.8), constrained_layout=True)

    ax1.hist(cs["cpbb"], bins=60, range=(0, 0.3), color=HEADLINE, alpha=0.85, label="real dark genes")
    ax1.axvline(d["cpbb"]["shuffle_median"], color=REF, ls=":", lw=1.5,
               label=f"shuffle median {d['cpbb']['shuffle_median']:.3f}")
    ax1.axvline(d["cpbb"]["markov1_median"], color=BASELINE, ls="--", lw=1.5,
               label=f"Markov-1 median {d['cpbb']['markov1_median']:.3f}")
    ax1.axvline(d["cpbb"]["real_median"], color="black", lw=1.5,
               label=f"real median {d['cpbb']['real_median']:.3f}")
    ax1.set_xlabel("codon-position base bias (CPBB)")
    ax1.set_ylabel("genes")
    ax1.set_title(f"Real vs null CPBB (n={d['n_genes']:,})")
    ax1.legend(frameon=False, fontsize=7)

    aurocs = [d["coding_vs_noise_auroc"]["real_vs_shuffle"], d["coding_vs_noise_auroc"]["real_vs_markov1"]]
    ax2.bar(["real vs\nshuffle", "real vs\nMarkov-1"], aurocs, color=[HEADLINE, BASELINE])
    ax2.axhline(0.5, color=REF, ls=":", lw=1, label="chance")
    ax2.set_ylim(0, 1)
    ax2.set_ylabel("coding-vs-noise AUROC")
    ax2.set_title("Dark genes separate cleanly from nulls")
    ax2.legend(frameon=False, fontsize=8)
    for i, v in enumerate(aurocs):
        ax2.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=9)
    fig.suptitle("Layer 5 — coding-structure vs shuffled/Markov-1 nulls (P4-D5)", fontsize=11)
    save(fig, "layer5_coding_structure_summary.png")


def layer2_prostt5() -> None:
    # Committed numbers (P4-D6 / README §7), same 60-family / 2,698-protein
    # matched eval for all three models -- no raw file available locally to
    # re-derive from (ProstT5's .npz is gitignored, lives on Angshuman's Mac).
    models = ["ESM-2 L22\n(sequence)", "ProstT5\ncentered (structure)", "Genos-m L9\n(genome, 300-fam ref)"]
    means = [0.993, 0.960, 0.739]
    medians = [0.998, 0.972, 0.795]
    colors = [HEADLINE, THIRD, BASELINE]

    fig, ax = plt.subplots(figsize=(6.5, 4.0), constrained_layout=True)
    x = np.arange(len(models))
    w = 0.35
    ax.bar(x - w / 2, means, w, label="mean", color=colors, alpha=0.95)
    ax.bar(x + w / 2, medians, w, label="median", color=colors, alpha=0.55)
    ax.axhline(0.5, color=REF, ls=":", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(models, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("held-out-family AUROC")
    ax.set_title("Layer 2 — ProstT5 vs sequence and genomic arms\n(same 60 families, 2,698 proteins, P4-D6)", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    for i, (m, med) in enumerate(zip(means, medians)):
        ax.text(i - w / 2, m + 0.015, f"{m:.3f}", ha="center", fontsize=7)
        ax.text(i + w / 2, med + 0.015, f"{med:.3f}", ha="center", fontsize=7)
    save(fig, "layer2_prostt5_summary.png")


if __name__ == "__main__":
    layer4_convergence()
    layer5_coding_structure()
    layer2_prostt5()
