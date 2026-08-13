#!/usr/bin/env python
"""CLI: label panel genome proteins dark-at-T0 / characterised-at-T1(proxy),
per P1-D3's operational definition with the Pfam-37 net-new-family proxy
standing in for InterPro-latest (see pfam_diff.py docstring for why).

Processes genomes in resumable batches, appending to the output CSV after
each batch, instead of one single-shot run that only writes at the very
end -- found necessary after the environment restarted mid-run three
times in a row, each time losing the full multi-hour scan with nothing
to show for it. A restart now costs at most one batch's worth of time
(the already-labeled genomes on disk are skipped on the next run), not
the whole panel.

Within a batch, all its genomes' proteins are still scanned together in
one combined hmmsearch call per Pfam release (not per-genome) -- that
batching is what made this fast at all; see hmmscan.py.

Usage:
    uv run python scripts/label_panel_proteins.py --release R207
    uv run python scripts/label_panel_proteins.py --release R207 --limit-genomes 50
    uv run python scripts/label_panel_proteins.py --release R207 --batch-size 25
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

PFAM35_HMM = RAW_ROOT / "pfam_35.0" / "Pfam-A.hmm.gz"
PFAM35_CLANS = RAW_ROOT / "pfam_35.0" / "Pfam-A.clans.tsv.gz"
PFAM37_HMM = RAW_ROOT / "pfam_37.0" / "Pfam-A.hmm.gz"
PFAM37_CLANS = RAW_ROOT / "pfam_37.0" / "Pfam-A.clans.tsv.gz"

CSV_HEADER = ["genome_accession", "protein_id", "dark_at_t0", "characterised_t1_proxy", "positive_proxy"]


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
            n_dark += row["dark_at_t0"] == "True"
            n_char += row["characterised_t1_proxy"] == "True"
            n_pos += row["positive_proxy"] == "True"
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

    pfam35_hits = scan_against_pfam(sequences, PFAM35_HMM, cpus=cpus)
    pfam37_hits = scan_against_pfam(sequences, PFAM37_HMM, cpus=cpus)

    write_header = not out_path.exists()
    with out_path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(CSV_HEADER)
        for protein_id in pfam35_hits:
            label = label_protein(pfam35_hits[protein_id], pfam37_hits.get(protein_id, []), new_family_ids)
            genome_accession = genome_of_protein.get(protein_id, "unknown")
            writer.writerow(
                [genome_accession, protein_id, label["dark_at_t0"], label["characterised_t1_proxy"], label["positive_proxy"]]
            )

    return len(sequences)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="R207")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--limit-genomes", type=int, default=None, help="cap for a quick pilot run")
    args = parser.parse_args()

    proteins_dir = PROC_ROOT / f"gtdb_{args.release}" / "panel_proteins"
    faa_files = sorted(proteins_dir.glob("*_protein.faa"))
    if args.limit_genomes:
        faa_files = faa_files[: args.limit_genomes]
    if not faa_files:
        raise SystemExit(f"no protein files in {proteins_dir} -- run extract_panel_proteins.py first")

    out_path = PROC_ROOT / f"gtdb_{args.release}" / "panel_protein_labels.csv"
    done = _already_labeled_genomes(out_path)
    remaining = [f for f in faa_files if f.name.removesuffix("_protein.faa") not in done]

    print(f"{len(done)}/{len(faa_files)} genomes already labeled, {len(remaining)} remaining", flush=True)
    if not remaining:
        n_total, n_dark, n_char, n_pos = _summarize(out_path)
        print(f"{n_total} proteins, dark-at-T0 {n_dark} ({100*n_dark/n_total:.2f}%), positive-proxy {n_pos} ({100*n_pos/n_total:.3f}%)")
        return

    print(f"loading pfam family diff (35.0 -> 37.0)...", flush=True)
    new_family_ids = new_families_since(PFAM35_CLANS, PFAM37_CLANS)

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
    print(f"dark-at-T0: {n_dark} ({100*n_dark/n_total:.2f}%)")
    print(f"characterised-T1-proxy: {n_char} ({100*n_char/n_total:.2f}%)")
    print(f"positive-proxy (dark-at-T0 AND characterised-T1-proxy): {n_pos} ({100*n_pos/n_total:.3f}%)")

    log_experiment(
        title=f"snapshot differencing ({args.release} panel, {len(faa_files)} genomes)",
        config=f"pfam-35 GA dark-at-T0 / pfam-37 net-new-family characterised-T1-proxy, batch-size={args.batch_size}",
        result=(
            f"{n_total} proteins; dark-at-T0 {n_dark} ({100*n_dark/n_total:.2f}%); "
            f"characterised-T1-proxy {n_char} ({100*n_char/n_total:.2f}%); "
            f"positive-proxy {n_pos} ({100*n_pos/n_total:.3f}%)"
        ),
        next_step="report to Track 1 for go/no-go sign-off" if n_pos >= 50 else "below the 50-100 floor, widen interval or broaden references per D6",
    )


if __name__ == "__main__":
    main()
