#!/usr/bin/env python
"""CLI: label panel genome proteins dark-at-T0 / characterised-at-T1(proxy),
per P1-D3's operational definition with the Pfam-37 net-new-family proxy
standing in for InterPro-latest (see pfam_diff.py docstring for why).

Loads every requested genome's proteins into one combined sequence set and
scans each Pfam library (35.0, 37.0) exactly once across the whole panel
-- scanning genome-by-genome instead re-parses the ~300MB+ HMM files for
every genome, which took ~9 minutes for a single genome in testing (would
be ~77 hours for the full 502-genome panel).

Usage:
    uv run python scripts/label_panel_proteins.py --release R207
    uv run python scripts/label_panel_proteins.py --release R207 --limit-genomes 50
"""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

from pyhmmer.easel import Alphabet

from darkmatter.data.hmmscan import load_protein_sequences_multi, scan_against_pfam
from darkmatter.data.pfam_diff import label_protein, new_families_since

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"
RAW_ROOT = Path(__file__).resolve().parents[1] / "data" / "raw"

PFAM35_HMM = RAW_ROOT / "pfam_35.0" / "Pfam-A.hmm.gz"
PFAM35_CLANS = RAW_ROOT / "pfam_35.0" / "Pfam-A.clans.tsv.gz"
PFAM37_HMM = RAW_ROOT / "pfam_37.0" / "Pfam-A.hmm.gz"
PFAM37_CLANS = RAW_ROOT / "pfam_37.0" / "Pfam-A.clans.tsv.gz"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--limit-genomes", type=int, default=None, help="cap for a quick pilot run")
    args = parser.parse_args()

    proteins_dir = PROC_ROOT / f"gtdb_{args.release}" / "panel_proteins"
    faa_files = sorted(proteins_dir.glob("*_protein.faa"))
    if args.limit_genomes:
        faa_files = faa_files[: args.limit_genomes]
    if not faa_files:
        raise SystemExit(f"no protein files in {proteins_dir} -- run extract_panel_proteins.py first")

    print(f"loading pfam family diff (35.0 -> 37.0)...", flush=True)
    new_family_ids = new_families_since(PFAM35_CLANS, PFAM37_CLANS)
    print(f"{len(new_family_ids)} net-new pfam families since 35.0", flush=True)

    print(f"loading proteins from {len(faa_files)} genomes...", flush=True)
    t0 = time.time()
    alphabet = Alphabet.amino()
    sequences = load_protein_sequences_multi(faa_files, alphabet)
    genome_of_protein: dict[str, str] = {}
    for faa_path in faa_files:
        genome_accession = faa_path.name.removesuffix("_protein.faa")
        with open(faa_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    genome_of_protein[line[1:].split(None, 1)[0]] = genome_accession
    print(f"{len(sequences)} proteins loaded in {time.time() - t0:.0f}s", flush=True)

    print("scanning against Pfam-35 (dark-at-T0)...", flush=True)
    t0 = time.time()
    pfam35_hits = scan_against_pfam(sequences, PFAM35_HMM, cpus=args.cpus)
    print(f"  done in {time.time() - t0:.0f}s", flush=True)

    print("scanning against Pfam-37 (characterised-at-T1 proxy)...", flush=True)
    t0 = time.time()
    pfam37_hits = scan_against_pfam(sequences, PFAM37_HMM, cpus=args.cpus)
    print(f"  done in {time.time() - t0:.0f}s", flush=True)

    out_path = PROC_ROOT / f"gtdb_{args.release}" / "panel_protein_labels.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_dark = n_characterised_proxy = n_positive_proxy = n_total = 0

    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["genome_accession", "protein_id", "dark_at_t0", "characterised_t1_proxy", "positive_proxy"])

        for protein_id in pfam35_hits:
            label = label_protein(pfam35_hits[protein_id], pfam37_hits.get(protein_id, []), new_family_ids)
            genome_accession = genome_of_protein.get(protein_id, "unknown")
            writer.writerow(
                [genome_accession, protein_id, label["dark_at_t0"], label["characterised_t1_proxy"], label["positive_proxy"]]
            )
            n_total += 1
            n_dark += label["dark_at_t0"]
            n_characterised_proxy += label["characterised_t1_proxy"]
            n_positive_proxy += label["positive_proxy"]

    print()
    print(f"wrote {out_path}")
    print(f"{n_total} proteins across {len(faa_files)} genomes")
    print(f"dark-at-T0: {n_dark} ({100*n_dark/n_total:.2f}%)")
    print(f"characterised-T1-proxy: {n_characterised_proxy} ({100*n_characterised_proxy/n_total:.2f}%)")
    print(f"positive-proxy (dark-at-T0 AND characterised-T1-proxy): {n_positive_proxy} ({100*n_positive_proxy/n_total:.3f}%)")
    if n_positive_proxy < 100:
        print(f"NOTE: {n_positive_proxy} positive-proxy sequences so far, floor is 50-100 (blueprint D6 go/no-go)")


if __name__ == "__main__":
    main()
