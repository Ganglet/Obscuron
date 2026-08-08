#!/usr/bin/env python
"""CLI: stream a GTDB representative-genome protein FASTA archive straight into S3.

Usage:
    uv run python scripts/stream_gtdb_to_s3.py --release R232 --bucket darkmatter-gtdb-067620369122
"""

from __future__ import annotations

import argparse

from darkmatter.data.s3_stream import stream_to_s3_with_retries

BASE_URL = "https://data.gtdb.ecogenomic.org/releases"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", required=True, help="e.g. R232")
    parser.add_argument("--bucket", required=True)
    parser.add_argument("--profile", default="track2")
    args = parser.parse_args()

    n = args.release.lstrip("Rr")
    url = f"{BASE_URL}/release{n}/{n}.0/genomic_files_reps/gtdb_proteins_aa_reps_r{n}.tar.gz"
    key = f"gtdb/{args.release}/genomic_files_reps/gtdb_proteins_aa_reps_r{n}.tar.gz"
    stream_to_s3_with_retries(url, args.bucket, key, args.profile)


if __name__ == "__main__":
    main()
