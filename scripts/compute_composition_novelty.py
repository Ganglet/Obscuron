#!/usr/bin/env python
"""Third convergence axis (Phase 4): model-free compositional/statistical novelty.

For each dark query, the Mahalanobis distance of its amino-acid composition from
the characterised-at-T0 reference's composition distribution. This is a genuinely
independent axis -- it shares no input with the ESM-2 embedding (Layer 1) or the
genomic context (Layer 4), and it is the buildable, hardware-free stand-in for the
cloud-blocked Genos-m axis, grounded in blueprint Layer 5 (statistical
discrimination) / the SS4 computational-linguistics parallel.

Emits composition_novelty.csv (protein_id, composition_novelty) -- fed to
run_convergence.py via --extra-axis for the 3-way convergence.

Usage:
    uv run python scripts/compute_composition_novelty.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from darkmatter.convergence import aa_composition, mahalanobis_novelty

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "gtdb_R207"
SAMPLE = PROC / "embedding_sample.csv"
PROTEINS = PROC / "panel_proteins"
OUT = PROC / "composition_novelty.csv"


def _load_sequences(df: pd.DataFrame) -> dict[str, str]:
    wanted = set(df["protein_id"])
    seqs: dict[str, str] = {}
    for acc in df["genome_accession"].unique():
        p = PROTEINS / f"{acc}_protein.faa"
        if not p.exists():
            continue
        for rec in SeqIO.parse(p, "fasta"):
            if rec.id in wanted:
                seqs[rec.id] = str(rec.seq)
    return seqs


def main() -> None:
    sample = pd.read_csv(SAMPLE)
    ref = sample[sample["category"] == "characterised_at_t0"]
    qry = sample[sample["category"].isin(["dark_negative", "positive"])]
    print(f"reference (characterised): {len(ref)}   queries (dark): {len(qry)}", flush=True)

    ref_seqs = _load_sequences(ref)
    qry_seqs = _load_sequences(qry)
    print(f"loaded sequences: {len(ref_seqs)} reference, {len(qry_seqs)} query", flush=True)

    ref_comp = np.array([aa_composition(ref_seqs[i]) for i in ref["protein_id"] if i in ref_seqs])
    q_ids = [i for i in qry["protein_id"] if i in qry_seqs]
    q_comp = np.array([aa_composition(qry_seqs[i]) for i in q_ids])

    novelty = mahalanobis_novelty(q_comp, ref_comp)
    out = pd.DataFrame({"protein_id": q_ids, "composition_novelty": novelty})
    out.to_csv(OUT, index=False)
    print(f"composition novelty for {len(out)} queries "
          f"(min {novelty.min():.2f}, median {np.median(novelty):.2f}, max {novelty.max():.2f})", flush=True)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
