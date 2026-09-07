#!/usr/bin/env python
"""Phase-1 robustness nice-to-have (P1-D2 §1): a second, shorter and more
recent snapshot boundary alongside the committed main boundary.

  main boundary   : Pfam 35.0 (Nov 2021) -> 37.0 (Jun 2024)   ~2.5 years
  robustness gap  : Pfam 36.0 (Sep 2023) -> 37.0 (Jun 2024)   ~9 months

The question this answers is the one a reviewer asks about any temporal split:
is the positive population an artifact of the specific (wide) boundary, or is
it robust to a different, shorter, more recent boundary? Reported at the two
levels that match the two Phase-1 positive estimates:

  1. Family-population diff -- parallels the P1-D5 go/no-go proxy, which was
     itself net-new-family based. Net-new Pfam families per window; the
     positive-generating mechanism is "a dark protein now hits a family that
     did not exist at T0", so net-new-family count is that mechanism's supply.

  2. Per-protein recount on a seeded subsample of the panel -- parallels the
     P1-D12 definitive per-protein count (4,138 positives, full panel). At the
     36->37 boundary: dark = no Pfam-36 hit at GA; characterised = hits a
     Pfam-37 family net-new since 36; positive = dark AND characterised. The
     35->37 numbers for the *same* genomes are read from the existing
     panel_protein_labels.csv, so the boundary is the only thing that changes.

Only the 36->37 boundary needs fresh hmmsearch (Pfam-36 for dark, Pfam-37 for
characterised); Pfam-35/37 main-boundary labels are reused from disk.

Usage:
    uv run python scripts/robustness_boundary.py                 # family diff + 60-genome subsample
    uv run python scripts/robustness_boundary.py --n-genomes 100
    uv run python scripts/robustness_boundary.py --family-only    # skip the per-protein scan
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import time
from pathlib import Path

from pyhmmer.easel import Alphabet

from darkmatter.data.hmmscan import load_protein_sequences_multi, scan_against_pfam
from darkmatter.data.pfam_diff import _strip_version, load_family_ids, new_families_since

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed" / "gtdb_R207"

PFAM35_CLANS = RAW / "pfam_35.0" / "Pfam-A.clans.tsv.gz"
PFAM36_CLANS = RAW / "pfam_36.0" / "Pfam-A.clans.tsv.gz"
PFAM37_CLANS = RAW / "pfam_37.0" / "Pfam-A.clans.tsv.gz"
PFAM36_HMM = RAW / "pfam_36.0" / "Pfam-A.hmm.gz"
PFAM37_HMM = RAW / "pfam_37.0" / "Pfam-A.hmm.gz"

PANEL_CSV = PROC / "genome_panel.csv"
LABELS_CSV = PROC / "panel_protein_labels.csv"
PROTEINS_DIR = PROC / "panel_proteins"
OUT_JSON = PROC / "robustness_boundary.json"


def family_population_diff() -> dict:
    f35 = load_family_ids(PFAM35_CLANS)
    f36 = load_family_ids(PFAM36_CLANS)
    f37 = load_family_ids(PFAM37_CLANS)
    return {
        "pfam35_families": len(f35),
        "pfam36_families": len(f36),
        "pfam37_families": len(f37),
        "net_new_35_37_gross": len(f37 - f35),           # families the label logic counts as "new since 35"
        "net_new_35_37_count_delta": len(f37) - len(f35),  # the P1-D5 headline (+2,347), net of retirements
        "net_new_35_36_gross": len(f36 - f35),
        "net_new_36_37_gross": len(f37 - f36),             # the shorter, recent window's supply
        "retired_35_37": len(f35 - f37),
    }


def _subsample_genomes(n: int, seed: int) -> list[str]:
    with PANEL_CSV.open(newline="", encoding="utf-8") as f:
        accs = [row["accession"] for row in csv.DictReader(f)]
    accs = [a for a in accs if (PROTEINS_DIR / f"{a}_protein.faa").exists()]
    rng = random.Random(seed)
    return sorted(rng.sample(accs, min(n, len(accs))))


def _main_boundary_labels(subsample: set[str]) -> dict[str, dict]:
    """Read the committed 35->37 labels for the subsample genomes' proteins."""
    out: dict[str, dict] = {}
    with LABELS_CSV.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["genome_accession"] in subsample:
                out[row["protein_id"]] = {
                    "dark35": row["dark_at_t0"] == "True",
                    "pos35_37": row["positive_proxy"] == "True",
                }
    return out


def per_protein_recount(n_genomes: int, seed: int, cpus: int) -> dict:
    subsample = set(_subsample_genomes(n_genomes, seed))
    faa_files = [PROTEINS_DIR / f"{a}_protein.faa" for a in sorted(subsample)]
    print(f"subsample: {len(subsample)} genomes (seed {seed})", flush=True)

    new_since_36 = new_families_since(PFAM36_CLANS, PFAM37_CLANS)

    alphabet = Alphabet.amino()
    sequences = load_protein_sequences_multi(faa_files, alphabet)
    n_proteins = len(sequences)
    print(f"loaded {n_proteins} proteins; scanning vs Pfam-36 then Pfam-37...", flush=True)

    t0 = time.time()
    hits36 = scan_against_pfam(sequences, PFAM36_HMM, cpus=cpus)
    print(f"  Pfam-36 scan: {time.time()-t0:.0f}s", flush=True)
    t1 = time.time()
    hits37 = scan_against_pfam(sequences, PFAM37_HMM, cpus=cpus)
    print(f"  Pfam-37 scan: {time.time()-t1:.0f}s", flush=True)

    main = _main_boundary_labels(subsample)

    n = dark36 = pos36_37 = 0
    dark35 = pos35_37 = 0
    for pid in hits36:
        n += 1
        is_dark36 = len(hits36[pid]) == 0
        is_char = any(_strip_version(a) in new_since_36 for a in hits37.get(pid, []))
        dark36 += is_dark36
        pos36_37 += is_dark36 and is_char
        m = main.get(pid)
        if m is not None:
            dark35 += m["dark35"]
            pos35_37 += m["pos35_37"]

    return {
        "n_genomes": len(subsample),
        "n_proteins": n,
        "seed": seed,
        # 36->37 shorter/recent boundary (fresh scans)
        "dark_at_36": dark36,
        "dark_frac_36": dark36 / n,
        "positive_36_37": pos36_37,
        "positive_rate_36_37": pos36_37 / n,
        # 35->37 main boundary, SAME genomes (from committed labels)
        "dark_at_35": dark35,
        "dark_frac_35": dark35 / n,
        "positive_35_37": pos35_37,
        "positive_rate_35_37": pos35_37 / n,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-genomes", type=int, default=60)
    ap.add_argument("--seed", type=int, default=42, help="frozen Phase-1 panel seed (P1-D6)")
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--family-only", action="store_true")
    args = ap.parse_args()

    result = {"family_population": family_population_diff()}
    fp = result["family_population"]
    print("=== family-population diff (positive-generating supply) ===")
    print(f"Pfam families: 35={fp['pfam35_families']}  36={fp['pfam36_families']}  37={fp['pfam37_families']}")
    print(f"net-new 35->37 (main, gross): {fp['net_new_35_37_gross']}  (P1-D5 net count +{fp['net_new_35_37_count_delta']})")
    print(f"net-new 36->37 (shorter/recent, gross): {fp['net_new_36_37_gross']}")
    print(f"  -> shorter window retains {100*fp['net_new_36_37_gross']/fp['net_new_35_37_gross']:.0f}% of the main window's new-family supply")

    if not args.family_only:
        print("\n=== per-protein recount (seeded panel subsample) ===")
        pp = per_protein_recount(args.n_genomes, args.seed, args.cpus)
        result["per_protein"] = pp
        print(f"\n{pp['n_proteins']} proteins across {pp['n_genomes']} genomes")
        print(f"dark-fraction : 35-boundary {100*pp['dark_frac_35']:.2f}%   36-boundary {100*pp['dark_frac_36']:.2f}%")
        print(f"positive-rate : 35->37 {100*pp['positive_rate_35_37']:.3f}%   36->37 {100*pp['positive_rate_36_37']:.3f}%")
        print(f"positives     : 35->37 {pp['positive_35_37']}   36->37 {pp['positive_36_37']} (on the subsample)")

    OUT_JSON.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nwrote {OUT_JSON}")


if __name__ == "__main__":
    main()
