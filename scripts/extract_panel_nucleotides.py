#!/usr/bin/env python
"""CLI: extract the genome panel's nucleotide CDS FASTA directly from
GTDB's public mirror -- the modality Genos-m needs (single-nucleotide
tokenization), not the amino-acid proteins already in S3/local disk.

No AWS credentials required -- streams straight from the public GTDB
mirror, not through S3.

Usage:
    uv run python scripts/extract_panel_nucleotides.py --release R207
    uv run python scripts/extract_panel_nucleotides.py --release R207 --accessions-from data/processed/gtdb_R207/embedding_sample.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from darkmatter.data.panel_nucleotides import extract_panel_nucleotides_from_gtdb

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207")
    parser.add_argument(
        "--accessions-from",
        type=Path,
        default=None,
        help="CSV with a genome_accession column to restrict to (default: the full genome_panel.csv)",
    )
    args = parser.parse_args()

    if args.accessions_from:
        with args.accessions_from.open(newline="", encoding="utf-8") as f:
            accessions = {row["genome_accession"] for row in csv.DictReader(f)}
    else:
        panel_csv = PROC_ROOT / f"gtdb_{args.release}" / "genome_panel.csv"
        with panel_csv.open(newline="", encoding="utf-8") as f:
            accessions = {row["accession"] for row in csv.DictReader(f)}

    out_dir = PROC_ROOT / f"gtdb_{args.release}" / "panel_nucleotides"
    existing = {p.name.removesuffix("_protein.fna") for p in out_dir.glob("*_protein.fna")} if out_dir.exists() else set()
    missing = accessions - existing

    if not missing:
        print(f"all {len(accessions)} panel nucleotide files already extracted in {out_dir}")
        return

    print(f"extracting {len(missing)}/{len(accessions)} missing panel nucleotide sequences from GTDB's public mirror", flush=True)
    written = extract_panel_nucleotides_from_gtdb(args.release, missing, out_dir)
    print(f"wrote {len(written)} nucleotide FASTA files to {out_dir}")

    still_missing = accessions - (existing | {p.name.removesuffix("_protein.fna") for p in written})
    if still_missing:
        print(f"{len(still_missing)} panel accessions were never found in the archive")


if __name__ == "__main__":
    main()
