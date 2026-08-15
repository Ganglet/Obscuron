#!/usr/bin/env python
"""CLI: embed the Stage 3 sample (build_embedding_sample.py output).

Defaults to ESM-2 -- the leakage-clean model for the retrospective
headline per P1-D7 (its UniRef50 2021_04 training data predates every
T0->T1 characterisation event; Genos-m's GTDB R220 pretraining sits
inside the window). Genos-m is not runnable on this machine's RTX 4060
(OOMs even quantized, see docs/reproducibility.md) -- pass --model
genos-m on a machine that can load it.

Usage:
    uv run python scripts/embed_panel_sample.py --release R207
"""

from __future__ import annotations

import time
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from darkmatter.embeddings import load_embedder
from darkmatter.experiment_log import log_experiment

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207")
    parser.add_argument("--model", choices=["genos-m", "esm2"], default="esm2")
    parser.add_argument("--batch-size", type=int, default=None)
    args = parser.parse_args()

    out_dir = PROC_ROOT / f"gtdb_{args.release}"
    sample = pd.read_csv(out_dir / "embedding_sample.csv")

    proteins_dir = out_dir / "panel_proteins"
    seq_by_id: dict[str, str] = {}
    print(f"loading sequences for {len(sample)} sampled proteins across "
          f"{sample['genome_accession'].nunique()} genomes...", flush=True)
    for genome_accession, group in sample.groupby("genome_accession"):
        faa_path = proteins_dir / f"{genome_accession}_protein.faa"
        wanted = set(group["protein_id"])
        for record in SeqIO.parse(faa_path, "fasta"):
            if record.id in wanted:
                seq_by_id[record.id] = str(record.seq)

    missing = sample.loc[~sample["protein_id"].isin(seq_by_id), "protein_id"].tolist()
    if missing:
        print(f"WARNING: {len(missing)} protein IDs not found in extracted FASTA, dropping", flush=True)
        sample = sample[~sample["protein_id"].isin(missing)].reset_index(drop=True)

    seqs = [seq_by_id[pid] for pid in sample["protein_id"]]

    print(f"loading {args.model}...", flush=True)
    embedder = load_embedder(args.model)
    kwargs = {} if args.batch_size is None else {"batch_size": args.batch_size}

    print(f"embedding {len(seqs)} sequences...", flush=True)
    t0 = time.time()
    vectors = embedder.embed(seqs, **kwargs)
    elapsed = time.time() - t0
    print(f"embedded {len(seqs)} sequences in {elapsed:.0f}s -> shape {vectors.shape}", flush=True)

    npy_path = out_dir / f"{args.model}_panel_embeddings.npy"
    manifest_path = out_dir / f"{args.model}_panel_embeddings_manifest.csv"
    np.save(npy_path, vectors)
    sample.to_csv(manifest_path, index=False)

    print(f"wrote {npy_path}")
    print(f"wrote {manifest_path}")

    counts = sample["category"].value_counts().to_dict()
    log_experiment(
        title=f"stage 3 embedding sample ({args.model}, {args.release} panel)",
        config=f"{args.model}, {len(sample)} sequences: {counts}",
        result=f"embedded in {elapsed:.0f}s, shape {vectors.shape}",
        next_step="feed into the Layer 1 novelty scorer once Track 1 fixes EVT vs density-based",
    )


if __name__ == "__main__":
    main()
