#!/usr/bin/env python
"""Layer 5 driver (Phase-4): statistical discrimination of coding structure from noise.

For each dark-matter gene, compute codon-position base bias (CPBB, from the
nucleotide CDS) + amino-acid k-mer entropy, then test whether these separate the
real genes from their own shuffled and Markov-1 nulls (the control group). A high
real-vs-null AUROC is evidence the dark genes carry genuine coding structure and
are not spurious ORF calls (blueprint SS6 Layer 5 / SS10 non-coding-artifact risk).

Emits coding_structure.csv (per-gene CPBB + entropy) for optional use as a
quality filter on the convergent set.

Usage:
    uv run python scripts/run_statistical.py
    uv run python scripts/run_statistical.py --max-genes 5000   # quick pass
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from darkmatter.statistical import (
    codon_position_bias,
    kmer_entropy,
    markov1_null,
    shuffle_seq,
)

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "gtdb_R207"
SAMPLE = PROC / "embedding_sample.csv"
AA_DIR = PROC / "panel_proteins"
NT_DIR = PROC / "panel_nucleotides"
OUT_CSV = PROC / "coding_structure.csv"
OUT_JSON = PROC / "statistical_results.json"


def _load(df: pd.DataFrame, directory: Path, suffix: str) -> dict[str, str]:
    wanted = set(df["protein_id"])
    seqs: dict[str, str] = {}
    for acc in df["genome_accession"].unique():
        p = directory / f"{acc}{suffix}"
        if p.exists():
            for rec in SeqIO.parse(p, "fasta"):
                if rec.id in wanted:
                    seqs[rec.id] = str(rec.seq)
    return seqs


def _auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    """AUROC that pos > neg (tie-safe via average ranks)."""
    pos = pos[~np.isnan(pos)]; neg = neg[~np.isnan(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank().to_numpy()
    return (ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-genes", type=int, default=0, help="0 = all dark queries")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    sample = pd.read_csv(SAMPLE)
    qry = sample[sample["category"].isin(["dark_negative", "positive"])].reset_index(drop=True)
    if args.max_genes:
        qry = qry.sample(n=min(args.max_genes, len(qry)), random_state=args.seed).reset_index(drop=True)
    print(f"{len(qry)} dark queries", flush=True)

    nt = _load(qry, NT_DIR, "_protein.fna")
    aa = _load(qry, AA_DIR, "_protein.faa")
    print(f"loaded {len(nt)} nt + {len(aa)} aa sequences", flush=True)

    rows, real_cpbb, shuf_cpbb, mark_cpbb = [], [], [], []
    for pid in qry["protein_id"]:
        if pid not in nt or pid not in aa:
            continue
        c_real = codon_position_bias(nt[pid])
        c_shuf = codon_position_bias(shuffle_seq(nt[pid], rng))
        c_mark = codon_position_bias(markov1_null(nt[pid], rng))
        real_cpbb.append(c_real); shuf_cpbb.append(c_shuf); mark_cpbb.append(c_mark)
        rows.append({"protein_id": pid, "cpbb": c_real,
                     "aa_entropy_k3": kmer_entropy(aa[pid], 3)})

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    real_cpbb = np.array(real_cpbb); shuf_cpbb = np.array(shuf_cpbb); mark_cpbb = np.array(mark_cpbb)

    auroc_shuf = _auroc(real_cpbb, shuf_cpbb)
    auroc_mark = _auroc(real_cpbb, mark_cpbb)
    result = {
        "n_genes": int(len(out)),
        "seed": args.seed,
        "cpbb": {"real_median": float(np.nanmedian(real_cpbb)),
                  "shuffle_median": float(np.nanmedian(shuf_cpbb)),
                  "markov1_median": float(np.nanmedian(mark_cpbb))},
        "coding_vs_noise_auroc": {"real_vs_shuffle": auroc_shuf, "real_vs_markov1": auroc_mark},
        "aa_entropy_k3_median": float(np.nanmedian(out["aa_entropy_k3"])),
    }
    OUT_JSON.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("\n=== coding-structure (CPBB) — real dark genes vs their nulls ===")
    print(f"CPBB median : real {result['cpbb']['real_median']:.4f}  "
          f"shuffle {result['cpbb']['shuffle_median']:.4f}  markov1 {result['cpbb']['markov1_median']:.4f}")
    print(f"coding-vs-noise AUROC : real vs shuffle {auroc_shuf:.4f}   real vs markov1 {auroc_mark:.4f}")
    print(f"(AUROC >> 0.5 = the dark genes carry genuine codon structure, not noise)")
    print(f"aa 3-mer entropy median: {result['aa_entropy_k3_median']:.3f} bits")
    print(f"\nwrote {OUT_CSV}\nwrote {OUT_JSON}")


if __name__ == "__main__":
    main()
