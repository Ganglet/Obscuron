#!/usr/bin/env python
"""CLI: download a GTDB or Pfam release by tag into data/raw/<source>_<release>/.

Usage:
    uv run python scripts/fetch_snapshot.py --source gtdb --release R226 --metadata-only
    uv run python scripts/fetch_snapshot.py --source pfam --release 37.0 --metadata-only
"""

from __future__ import annotations

import argparse
from pathlib import Path

from darkmatter.data.gtdb import fetch_gtdb_release
from darkmatter.data.manifest import write_manifest
from darkmatter.data.pfam import fetch_pfam_release

DATA_ROOT = Path(__file__).resolve().parents[1] / "data" / "raw"
MANIFEST_PATH = Path(__file__).resolve().parents[1] / "data" / "manifest.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=["gtdb", "pfam"], required=True)
    parser.add_argument("--release", required=True, help="e.g. R226 for gtdb, 37.0 for pfam")
    parser.add_argument("--metadata-only", action="store_true", help="fetch small metadata/taxonomy files only")
    args = parser.parse_args()

    if args.source == "gtdb":
        files = fetch_gtdb_release(args.release, DATA_ROOT, metadata_only=args.metadata_only)
    else:
        files = fetch_pfam_release(args.release, DATA_ROOT, metadata_only=args.metadata_only)

    write_manifest(MANIFEST_PATH, args.source, args.release, files, extra={"metadata_only": args.metadata_only})

    print(f"Fetched {len(files)} file(s) for {args.source} {args.release}:")
    for f in files:
        size_mb = f.stat().st_size / (1024**2)
        print(f"  {f}  ({size_mb:.2f} MB)")


if __name__ == "__main__":
    main()
