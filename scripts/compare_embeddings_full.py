#!/usr/bin/env python
"""CLI: embedding-separation comparison across a phylum-diverse sample of
panel genomes, not just the single-genome pilot (compare_embeddings_pilot.py).

Genos-m OOMs on this machine's RTX 4060 (see docs/reproducibility.md);
this runs ESM-2 only. Genos-m needs to run on Angshuman's M1 Pro, or this
same approach pointed at Genos-m there.

Usage:
    uv run python scripts/compare_embeddings_full.py --n-genomes 30
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from pyhmmer.easel import Alphabet

from darkmatter.data.hmmscan import load_protein_sequences_multi, scan_against_pfam
from darkmatter.embeddings import load_embedder
from darkmatter.experiment_log import log_experiment
from darkmatter.separation import select_labeled_subset, separation_score

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"
RAW_ROOT = Path(__file__).resolve().parents[1] / "data" / "raw"
PFAM35_HMM = RAW_ROOT / "pfam_35.0" / "Pfam-A.hmm.gz"


def _sample_genomes(release: str, n_genomes: int, seed: int) -> list[str]:
    """Phylum-stratified sample of accessions already present in genome_panel.csv,
    so the comparison spans real taxonomic diversity, not just whatever sorts first."""
    panel = pd.read_csv(PROC_ROOT / f"gtdb_{release}" / "genome_panel.csv")
    rng = np.random.default_rng(seed)
    per_phylum = max(1, n_genomes // panel["phylum"].nunique())
    chosen = (
        panel.groupby("phylum", group_keys=False)
        .apply(lambda g: g.sample(n=min(per_phylum, len(g)), random_state=rng))
    )
    if len(chosen) > n_genomes:
        chosen = chosen.sample(n=n_genomes, random_state=rng)
    return chosen["accession"].tolist()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207")
    parser.add_argument("--model", choices=["genos-m", "esm2"], default="esm2")
    parser.add_argument("--n-genomes", type=int, default=30)
    parser.add_argument("--per-family", type=int, default=5)
    parser.add_argument("--max-families", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cpus", type=int, default=4)
    args = parser.parse_args()

    proteins_dir = PROC_ROOT / f"gtdb_{args.release}" / "panel_proteins"
    accessions = _sample_genomes(args.release, args.n_genomes, args.seed)
    faa_files = [proteins_dir / f"{acc}_protein.faa" for acc in accessions]
    faa_files = [f for f in faa_files if f.exists()]
    print(f"sampled {len(faa_files)} genomes across phyla for the labeled subset", flush=True)

    alphabet = Alphabet.amino()
    sequences = load_protein_sequences_multi(faa_files, alphabet)
    print(f"{len(sequences)} proteins loaded", flush=True)

    print("scanning against Pfam-35 for family labels...", flush=True)
    t0 = time.time()
    hits = scan_against_pfam(sequences, PFAM35_HMM, cpus=args.cpus)
    print(f"  done in {time.time() - t0:.0f}s", flush=True)

    families = select_labeled_subset(hits, per_family=args.per_family, max_families=args.max_families)
    total = sum(len(v) for v in families.values())
    print(f"labeled subset: {len(families)} families, {total} sequences", flush=True)

    seq_by_id: dict[str, str] = {}
    for faa_path in faa_files:
        for record in SeqIO.parse(faa_path, "fasta"):
            seq_by_id[record.id] = str(record.seq)

    ordered_ids = [sid for ids in families.values() for sid in ids]
    labels = [fam for fam, ids in families.items() for _ in ids]
    seqs = [seq_by_id[sid] for sid in ordered_ids]

    print(f"loading {args.model}...", flush=True)
    embedder = load_embedder(args.model)
    t0 = time.time()
    vectors = embedder.embed(seqs)
    print(f"  embedded {len(seqs)} sequences in {time.time() - t0:.0f}s", flush=True)

    result = separation_score(vectors, labels)
    result["model"] = args.model
    result["n_genomes_sampled"] = len(faa_files)
    result["families"] = {fam: ids for fam, ids in families.items()}

    out_path = PROC_ROOT / f"gtdb_{args.release}" / f"{args.model}_separation_full.json"
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "families"}, indent=2))
    print(f"wrote {out_path}")

    log_experiment(
        title=f"{args.model} embedding separation ({len(faa_files)} genomes, {len(families)} families)",
        config=f"{args.model}, mean-pooled, unambiguous single-Pfam-35-hit proteins, {args.per_family}/family, seed={args.seed}",
        result=(
            f"{result['n_sequences']} sequences, {result['n_families']} families: "
            f"within-family cosine {result['within_family_mean_cosine']:.3f} vs "
            f"across-family {result['across_family_mean_cosine']:.3f}, "
            f"gap {result['separation_gap']:.3f}"
        ),
        next_step="compare against the other model's separation once both are available",
    )


if __name__ == "__main__":
    main()
