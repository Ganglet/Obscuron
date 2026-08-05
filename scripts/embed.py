#!/usr/bin/env python
"""CLI: embed sequences from a FASTA file with genos-m or esm2.

This is the concrete Track 2 deliverable for blueprint §8 Phase 1:
"Run the Genos-m/ESM-2 comparison on a small labeled subset and report
the result to Track 1." Output is a .npy embedding matrix plus a
sidecar .ids.txt with one sequence ID per row, in order.

Usage:
    uv run python scripts/embed.py --model genos-m --fasta path/to.fasta --out path/to/genos_m.npy
    uv run python scripts/embed.py --model esm2 --fasta path/to.fasta --out path/to/esm2.npy
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
from Bio import SeqIO

from darkmatter.embeddings import load_embedder


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=["genos-m", "esm2"], required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=None)
    args = parser.parse_args()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    if not records:
        raise SystemExit(f"No sequences found in {args.fasta}")
    ids = [r.id for r in records]
    seqs = [str(r.seq) for r in records]

    print(f"Loading {args.model}...")
    t0 = time.time()
    embedder = load_embedder(args.model)
    print(f"  loaded in {time.time() - t0:.1f}s, embedding_dim={embedder.embedding_dim}")

    kwargs = {} if args.batch_size is None else {"batch_size": args.batch_size}
    t0 = time.time()
    vectors = embedder.embed(seqs, **kwargs)
    print(f"Embedded {len(seqs)} sequences in {time.time() - t0:.1f}s -> shape {vectors.shape}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.out, vectors)
    args.out.with_suffix(".ids.txt").write_text("\n".join(ids), encoding="utf-8")
    print(f"Wrote {args.out} and {args.out.with_suffix('.ids.txt')}")


if __name__ == "__main__":
    main()
