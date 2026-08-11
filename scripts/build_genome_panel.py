#!/usr/bin/env python
"""CLI: build the phylum-stratified genome panel for the Phase 1 retrospective
benchmark (P1-D6, docs/Track1_phase1_benchmark_scope.md section 6).

Samples GTDB species-representative genomes from the historical (T0)
release, stratified by phylum with a floor/cap so rare phyla always
appear and dominant ones don't crowd out the panel.

Usage:
    uv run python scripts/build_genome_panel.py --release R207
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from darkmatter.data.panel import load_representative_accessions, sample_panel
from darkmatter.data.preprocess import load_taxonomy_release

DATA_ROOT = Path(__file__).resolve().parents[1] / "data" / "raw"
PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207", help="historical (T0) release the panel is drawn from")
    parser.add_argument("--target-total", type=int, default=500)
    parser.add_argument("--floor", type=int, default=2)
    parser.add_argument("--cap", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    taxonomy = load_taxonomy_release(args.release, DATA_ROOT)
    reps = load_representative_accessions(args.release, DATA_ROOT)
    panel, info = sample_panel(
        taxonomy,
        reps,
        target_total=args.target_total,
        floor=args.floor,
        cap=args.cap,
        seed=args.seed,
    )

    out_dir = PROC_ROOT / f"gtdb_{args.release}"
    out_dir.mkdir(parents=True, exist_ok=True)
    panel_path = out_dir / "genome_panel.csv"
    summary_path = out_dir / "genome_panel_summary.json"
    panel.to_csv(panel_path, index=False)
    summary_path.write_text(json.dumps(info, indent=2), encoding="utf-8")

    print(f"{args.release}: {len(reps)} species-rep genomes across {info['phyla_available']} phyla")
    print(f"panel: {info['panel_size']} genomes (target {args.target_total}), allocation f={info['allocation_f']:.4f}")
    print(f"phyla in panel: {info['phyla_in_panel']} / {info['phyla_available']}")
    print(f"wrote {panel_path}")
    print(f"wrote {summary_path}")
    print()
    print("NOTE: panel rows are accession + taxonomy only. Per-genome protein")
    print("FASTA source URLs (needed for Stage 2 hmmsearch labeling) aren't")
    print("resolvable from GTDB's taxonomy/metadata files alone — neither has")
    print("a download-path column. Needs a join against NCBI's assembly_summary")
    print("(or a decision from Track 1 on the intended source) before Stage 2.")


if __name__ == "__main__":
    main()
