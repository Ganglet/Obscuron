#!/usr/bin/env python
"""CLI: label panel genome proteins dark-at-Pfam36 / characterised-since-Pfam36,
at FULL panel scale (all 502 genomes) -- the Layer-4 Genos-m follow-up (see
docs/problems_and_decisions.md P4-D11).

Why this exists: the Layer-4 Genos-m convergence axis enriches for positives
(lift 1.30x, P4-D8) instead of depleting like every other axis, plausibly
because Genos-m's pretraining (GTDB R220, Apr 2024) sits inside the main
35.0->37.0 (Nov 2021 -> Jun 2024) boundary window -- genes characterised
anywhere in that ~2.5-year window could be leakage-contaminated. The tightest
boundary this project's already-fetched Pfam releases can support is
36.0 (Sep 2023) -> 37.0 (Jun 2024), ~9 months -- P1-D14's "robustness"
boundary, but that was only ever run on a 60-genome subsample. This script
runs the SAME 36->37 boundary at FULL PANEL SCALE (502 genomes) so the
leakage-reduced positive set used to re-test the Genos-m axis is a real
count, not an extrapolation from a subsample.

NOTE ON WHAT "LEAKAGE-CLEAN" MEANS HERE: Pfam 37.0 (~Jun 2024) sits only
~2 months after Genos-m's R220 cutoff (Apr 2024) -- there is no later Pfam
release fetched (or, as far as this project has confirmed, available) to
push T1 further out. So this is NOT a true "T0 >= R220" split; it is the
TIGHTEST available approximation given real Pfam release granularity. This
is a documented scoping limitation (same pattern as P1-D6/P2-D10/P3-D7/P4-D9),
not an oversight -- see the write-up in problems_and_decisions.md.

Processes genomes in resumable batches (same pattern as label_panel_proteins.py,
built after that script lost full multi-hour runs to session restarts three
times) -- a restart costs at most one batch, not the whole panel.

Usage:
    uv run python scripts/label_panel_proteins_36_37.py
    uv run python scripts/label_panel_proteins_36_37.py --batch-size 25 --cpus 4
"""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

from pyhmmer.easel import Alphabet

from darkmatter.data.hmmscan import load_protein_sequences_multi, scan_against_pfam
from darkmatter.data.pfam_diff import label_protein, new_families_since
from darkmatter.experiment_log import log_experiment

PROC_ROOT = Path(__file__).resolve().parents[1] / "data" / "processed"
RAW_ROOT = Path(__file__).resolve().parents[1] / "data" / "raw"

PFAM36_HMM = RAW_ROOT / "pfam_36.0" / "Pfam-A.hmm.gz"
PFAM36_CLANS = RAW_ROOT / "pfam_36.0" / "Pfam-A.clans.tsv.gz"
PFAM37_HMM = RAW_ROOT / "pfam_37.0" / "Pfam-A.hmm.gz"
PFAM37_CLANS = RAW_ROOT / "pfam_37.0" / "Pfam-A.clans.tsv.gz"

CSV_HEADER = ["genome_accession", "protein_id", "dark_at_36", "characterised_37_since_36", "positive_36_37"]


def _already_labeled_genomes(out_path: Path) -> set[str]:
    if not out_path.exists():
        return set()
    with out_path.open(newline="", encoding="utf-8") as f:
        return {row["genome_accession"] for row in csv.DictReader(f)}


def _summarize(out_path: Path) -> tuple[int, int, int, int]:
    n_total = n_dark = n_char = n_pos = 0
    with out_path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            n_total += 1
            n_dark += row["dark_at_36"] == "True"
            n_char += row["characterised_37_since_36"] == "True"
            n_pos += row["positive_36_37"] == "True"
    return n_total, n_dark, n_char, n_pos


def _label_batch(faa_files: list[Path], new_family_ids: set[str], cpus: int, out_path: Path) -> int:
    alphabet = Alphabet.amino()
    sequences = load_protein_sequences_multi(faa_files, alphabet)
    genome_of_protein: dict[str, str] = {}
    for faa_path in faa_files:
        genome_accession = faa_path.name.removesuffix("_protein.faa")
        with open(faa_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    genome_of_protein[line[1:].split(None, 1)[0]] = genome_accession

    pfam36_hits = scan_against_pfam(sequences, PFAM36_HMM, cpus=cpus)
    pfam37_hits = scan_against_pfam(sequences, PFAM37_HMM, cpus=cpus)

    write_header = not out_path.exists()
    with out_path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(CSV_HEADER)
        for protein_id in pfam36_hits:
            label = label_protein(pfam36_hits[protein_id], pfam37_hits.get(protein_id, []), new_family_ids)
            genome_accession = genome_of_protein.get(protein_id, "unknown")
            writer.writerow(
                [genome_accession, protein_id, label["dark_at_t0"], label["characterised_t1_proxy"], label["positive_proxy"]]
            )

    return len(sequences)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--limit-genomes", type=int, default=None, help="cap for a quick pilot run")
    args = parser.parse_args()

    proteins_dir = PROC_ROOT / "gtdb_R207" / "panel_proteins"
    faa_files = sorted(proteins_dir.glob("*_protein.faa"))
    if args.limit_genomes:
        faa_files = faa_files[: args.limit_genomes]
    if not faa_files:
        raise SystemExit(f"no protein files in {proteins_dir}")

    out_path = PROC_ROOT / "gtdb_R207" / "panel_protein_labels_36_37.csv"
    done = _already_labeled_genomes(out_path)
    remaining = [f for f in faa_files if f.name.removesuffix("_protein.faa") not in done]

    print(f"{len(done)}/{len(faa_files)} genomes already labeled, {len(remaining)} remaining", flush=True)
    if not remaining:
        n_total, n_dark, n_char, n_pos = _summarize(out_path)
        print(f"{n_total} proteins, dark-at-36 {n_dark} ({100*n_dark/n_total:.2f}%), positive-36-37 {n_pos} ({100*n_pos/n_total:.3f}%)")
        return

    print("loading pfam family diff (36.0 -> 37.0)...", flush=True)
    new_family_ids = new_families_since(PFAM36_CLANS, PFAM37_CLANS)

    n_batches = (len(remaining) + args.batch_size - 1) // args.batch_size
    for i in range(0, len(remaining), args.batch_size):
        batch = remaining[i : i + args.batch_size]
        batch_num = i // args.batch_size + 1
        t0 = time.time()
        n_proteins = _label_batch(batch, new_family_ids, args.cpus, out_path)
        elapsed = time.time() - t0
        print(f"batch {batch_num}/{n_batches}: {len(batch)} genomes, {n_proteins} proteins, {elapsed:.0f}s", flush=True)

    n_total, n_dark, n_char, n_pos = _summarize(out_path)
    print()
    print(f"wrote {out_path}")
    print(f"{n_total} proteins across {len(faa_files)} genomes")
    print(f"dark-at-36: {n_dark} ({100*n_dark/n_total:.2f}%)")
    print(f"characterised-37-since-36: {n_char} ({100*n_char/n_total:.2f}%)")
    print(f"positive-36-37 (dark-at-36 AND characterised-since-36): {n_pos} ({100*n_pos/n_total:.3f}%)")

    log_experiment(
        title="full-panel 36->37 boundary labeling (leakage-tightened, Genos-m axis follow-up)",
        config=f"pfam-36 GA dark-at-36 / pfam-37 net-new-family characterised-since-36, batch-size={args.batch_size}, full 502-genome panel",
        result=(
            f"{n_total} proteins; dark-at-36 {n_dark} ({100*n_dark/n_total:.2f}%); "
            f"characterised-37-since-36 {n_char} ({100*n_char/n_total:.2f}%); "
            f"positive-36-37 {n_pos} ({100*n_pos/n_total:.3f}%)"
        ),
        next_step="cross-reference against genos-m_novelty_scores.csv to recompute the convergence-axis lift on this leakage-tightened positive set",
    )


if __name__ == "__main__":
    main()
