#!/usr/bin/env python
"""CLI: parse a fetched GTDB release's taxonomy into a structured table + summary stats.

Usage:
    uv run python scripts/preprocess_taxonomy.py --release R232
    uv run python scripts/preprocess_taxonomy.py --release R207
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from darkmatter.data.preprocess import load_taxonomy_release, summarize

DATA_ROOT = Path(__file__).resolve().parents[1] / "data" / "raw"
OUT_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", required=True, help="e.g. R232")
    args = parser.parse_args()

    df = load_taxonomy_release(args.release, DATA_ROOT)
    stats = summarize(df)

    out_dir = OUT_ROOT / f"gtdb_{args.release}"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "taxonomy.csv", index=False)
    (out_dir / "summary.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")

    print(f"{args.release}:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
    print(f"wrote {out_dir / 'taxonomy.csv'}")


if __name__ == "__main__":
    main()
