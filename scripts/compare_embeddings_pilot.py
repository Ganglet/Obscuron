#!/usr/bin/env python
"""CLI: pilot embedding-separation check on one panel genome's real,
unambiguous Pfam-35 hits — the small-labeled-subset comparison from
blueprint section 8, run at pilot scale (one genome) ahead of the full
panel.

Genos-m OOMs on this machine's RTX 4060 (see docs/reproducibility.md);
this runs ESM-2 only. Genos-m needs to run on Angshuman's M1 Pro, or
this same script pointed at Genos-m there.

Usage:
    uv run python scripts/compare_embeddings_pilot.py \\
        --faa data/processed/gtdb_R207/panel_proteins/GB_GCA_001515205.2_protein.faa
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio import SeqIO

from darkmatter.data.hmmscan import scan_genome_against_pfam
from darkmatter.embeddings import load_embedder
from darkmatter.separation import select_labeled_subset, separation_score

PFAM_HMM = Path("data/raw/pfam_35.0/Pfam-A.hmm.gz")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--faa", type=Path, required=True)
    parser.add_argument("--model", choices=["genos-m", "esm2"], default="esm2")
    parser.add_argument("--per-family", type=int, default=4)
    parser.add_argument("--max-families", type=int, default=6)
    args = parser.parse_args()

    print(f"scanning {args.faa} against Pfam-35...", flush=True)
    hits = scan_genome_against_pfam(args.faa, PFAM_HMM, cpus=4)

    families = select_labeled_subset(hits, per_family=args.per_family, max_families=args.max_families)
    total = sum(len(v) for v in families.values())
    print(f"labeled subset: {len(families)} families, {total} sequences")
    for fam, ids in families.items():
        print(f"  {fam}: {ids}")

    seq_by_id = {r.id: str(r.seq) for r in SeqIO.parse(args.faa, "fasta")}
    ordered_ids = [sid for ids in families.values() for sid in ids]
    labels = [fam for fam, ids in families.items() for _ in ids]
    seqs = [seq_by_id[sid] for sid in ordered_ids]

    print(f"loading {args.model}...", flush=True)
    embedder = load_embedder(args.model)
    vectors = embedder.embed(seqs)

    result = separation_score(vectors, labels)
    result["model"] = args.model
    result["source_faa"] = str(args.faa)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
