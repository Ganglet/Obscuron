#!/usr/bin/env python
"""CLI: assemble the "protein table" Track 1's scoring module needs
(docs/Track2_Phase2_scoring_handoff.md) -- joins the embedding manifest
with GTDB phylum, Pfam-35 family labels, and protein length.

pfam_family is only resolved for reference (characterised-at-T0) rows --
a fresh Pfam-35 scan, since the Stage-2 labeling only kept a boolean
dark/characterised flag, not which family. Query (dark-at-T0) rows have
no Pfam-35 hit by definition, so pfam_family stays null for them; only
unambiguous (single-hit) reference proteins get a family id, since a
multi-family protein can't unambiguously be "the" held-out unit.

Usage:
    uv run python scripts/build_protein_table.py --release R207 --model esm2
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from pyhmmer.easel import Alphabet, DigitalSequenceBlock

from darkmatter.data.hmmscan import load_protein_sequences_multi, scan_against_pfam

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"
RAW_ROOT = Path(__file__).resolve().parents[1] / "data" / "raw"
PFAM35_HMM = RAW_ROOT / "pfam_35.0" / "Pfam-A.hmm.gz"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207")
    parser.add_argument("--model", choices=["genos-m", "esm2"], default="esm2")
    parser.add_argument("--cpus", type=int, default=4)
    args = parser.parse_args()

    out_dir = PROC_ROOT / f"gtdb_{args.release}"
    manifest = pd.read_csv(out_dir / f"{args.model}_panel_embeddings_manifest.csv")
    manifest["emb_row"] = range(len(manifest))
    manifest["role"] = manifest["category"].map(
        {"characterised_at_t0": "reference", "dark_negative": "query", "positive": "query"}
    )
    manifest["is_positive"] = manifest["category"] == "positive"

    panel = pd.read_csv(out_dir / "genome_panel.csv")
    manifest = manifest.merge(
        panel[["accession", "phylum"]], left_on="genome_accession", right_on="accession", how="left"
    ).drop(columns=["accession"])

    proteins_dir = out_dir / "panel_proteins"
    lengths: dict[str, int] = {}
    print(f"reading sequence lengths for {len(manifest)} proteins...", flush=True)
    for genome_accession, group in manifest.groupby("genome_accession"):
        faa_path = proteins_dir / f"{genome_accession}_protein.faa"
        if not faa_path.exists():
            continue
        wanted = set(group["protein_id"])
        for record in SeqIO.parse(faa_path, "fasta"):
            if record.id in wanted:
                lengths[record.id] = len(record.seq)
    manifest["length_aa"] = manifest["protein_id"].map(lengths)
    missing_length = manifest["length_aa"].isna().sum()
    if missing_length:
        print(f"WARNING: {missing_length} proteins missing a length (FASTA not found)", flush=True)

    ref_mask = manifest["role"] == "reference"
    ref_ids = set(manifest.loc[ref_mask, "protein_id"])
    ref_genomes = set(manifest.loc[ref_mask, "genome_accession"])
    faa_files = [proteins_dir / f"{acc}_protein.faa" for acc in ref_genomes]
    faa_files = [f for f in faa_files if f.exists()]

    alphabet = Alphabet.amino()
    all_sequences = load_protein_sequences_multi(faa_files, alphabet)
    ref_sequences = [s for s in all_sequences if s.name in ref_ids]
    ref_block = DigitalSequenceBlock(alphabet, ref_sequences)

    print(f"scanning {len(ref_block)} reference proteins against Pfam-35 for family labels...", flush=True)
    hits = scan_against_pfam(ref_block, PFAM35_HMM, cpus=args.cpus)
    pfam_family = {pid: (fams[0] if len(fams) == 1 else None) for pid, fams in hits.items()}
    manifest["pfam_family"] = manifest["protein_id"].map(pfam_family)

    out_path = out_dir / f"{args.model}_protein_table.csv"
    manifest.to_csv(out_path, index=False)

    n_with_family = manifest["pfam_family"].notna().sum()
    print(f"wrote {out_path}", flush=True)
    print(
        f"{len(manifest)} rows: {ref_mask.sum()} reference ({n_with_family} with unambiguous pfam_family), "
        f"{(~ref_mask).sum()} query",
        flush=True,
    )


if __name__ == "__main__":
    main()
