#!/usr/bin/env python
"""CLI: extract the genome panel's protein FASTA from the S3-hosted combined
GTDB rep-protein tarball, without downloading the full archive.

Usage:
    uv run python scripts/extract_panel_proteins.py --release R207 \\
        --bucket darkmatter-gtdb-067620369122
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from darkmatter.data.panel_proteins import extract_panel_proteins_from_s3

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207", help="panel's source release")
    parser.add_argument("--bucket", required=True)
    parser.add_argument("--profile", default="track2")
    args = parser.parse_args()

    n = args.release.lstrip("Rr")
    panel_csv = PROC_ROOT / f"gtdb_{args.release}" / "genome_panel.csv"
    with panel_csv.open(newline="", encoding="utf-8") as f:
        accessions = {row["accession"] for row in csv.DictReader(f)}

    out_dir = PROC_ROOT / f"gtdb_{args.release}" / "panel_proteins"
    existing = {p.name.removesuffix("_protein.faa") for p in out_dir.glob("*_protein.faa")} if out_dir.exists() else set()
    missing = accessions - existing

    if not missing:
        print(f"all {len(accessions)} panel proteins already extracted in {out_dir}")
        return

    key = f"gtdb/{args.release}/genomic_files_reps/gtdb_proteins_aa_reps_r{n}.tar.gz"
    print(f"extracting {len(missing)}/{len(accessions)} missing panel proteins from s3://{args.bucket}/{key}", flush=True)

    written = extract_panel_proteins_from_s3(args.bucket, key, args.profile, missing, out_dir)
    print(f"wrote {len(written)} protein FASTA files to {out_dir}")

    still_missing = accessions - (existing | {p.name.removesuffix("_protein.faa") for p in written})
    if still_missing:
        print(f"{len(still_missing)} panel accessions were never found in the archive")


if __name__ == "__main__":
    main()
