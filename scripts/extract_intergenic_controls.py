#!/usr/bin/env python
"""True-intergenic negative control for Layer 5 (coding-structure statistics).

Shuffled-sequence and Markov-1 nulls (src/darkmatter/statistical/) are
model-free substitutes for a real non-coding control; this script builds
the real thing -- actual intergenic DNA from the panel genomes -- by
streaming GTDB's full genome-assembly archive and extracting the sequence
between consecutive gene calls, using coordinates already embedded in the
Prodigal-style FASTA headers this project already has locally
(`>{contig}_{n} # {start} # {end} # {strand} # ...`), so no new
coordinate-calling step is needed.

Scope decision, not a shortcut (matches P1-D6/P2-D10's pattern): the
archive is 65GB and organised by accession-number path, not front-loaded
by panel membership, so covering all 502 panel genomes would mean
streaming close to the full archive. This streams a bounded prefix
(time- and byte-capped) and uses whichever panel genomes it encounters --
enough for a robust intergenic sample (thousands of gaps per genome) well
short of a full-archive read. Streamed directly to local processing, never
warehoused (P1-D4): only the small derived CPBB statistics are persisted,
not the genome sequences themselves.

    uv run python scripts/extract_intergenic_controls.py
"""
from __future__ import annotations

import json
import re
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "gtdb_R207"
NT_DIR = PROC / "panel_nucleotides"
OUT = PROC / "intergenic_controls.csv"
RESULTS_OUT = PROC / "statistical_results_intergenic.json"

ARCHIVE_URL = "https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/gtdb_genomes_reps_r207.tar.gz"
import os

MIN_GAP = 60          # shortest intergenic segment worth scoring
MAX_GENOMES = int(os.environ.get("INTERGENIC_MAX_GENOMES", 80))
MAX_BYTES = int(os.environ.get("INTERGENIC_MAX_BYTES", 12_000_000_000))     # ~12GB compressed-stream budget
MAX_SECONDS = int(os.environ.get("INTERGENIC_MAX_SECONDS", 1800))          # 30 min wall-clock budget
HEADER_RE = re.compile(r"^>(\S+)_(\d+) # (\d+) # (\d+) # ")


def panel_accession_map() -> dict[str, str]:
    """raw NCBI accession (as used in the archive path) -> panel accession
    (as used in our local filenames, e.g. GB_GCA_... / RS_GCF_...)."""
    panel = pd.read_csv(PROC / "genome_panel.csv")
    out = {}
    for acc in panel["accession"]:
        raw = acc.split("_", 1)[1] if acc.startswith(("GB_", "RS_")) else acc
        out[raw] = acc
    return out


def gene_coords(panel_accession: str) -> dict[str, list[tuple[int, int]]]:
    """contig -> sorted [(start, end), ...] parsed from the already-local
    nucleotide FASTA headers (Prodigal coordinates, 1-based inclusive)."""
    path = NT_DIR / f"{panel_accession}_protein.fna"
    coords: dict[str, list[tuple[int, int]]] = {}
    if not path.exists():
        return coords
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            m = HEADER_RE.match(line)
            if not m:
                continue
            contig, _n, start, end = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
            coords.setdefault(contig, []).append((start, end))
    for c in coords:
        coords[c].sort()
    return coords


def parse_fasta(text: str) -> dict[str, str]:
    seqs, cur, parts = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(parts)
            cur = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(parts)
    return seqs


def intergenic_segments(contig_seq: str, genes: list[tuple[int, int]]) -> list[str]:
    out = []
    prev_end = 0
    for start, end in genes:
        if start - 1 - prev_end >= MIN_GAP:
            out.append(contig_seq[prev_end:start - 1])
        prev_end = max(prev_end, end)
    if len(contig_seq) - prev_end >= MIN_GAP:
        out.append(contig_seq[prev_end:])
    return out


class _CountingRaw:
    """Wraps a streamed response's raw file object, tracking bytes read so
    we can enforce a compressed-byte budget mid-stream."""
    def __init__(self, raw):
        self._raw = raw
        self.total = 0

    def read(self, amt=-1):
        chunk = self._raw.read(amt)
        self.total += len(chunk)
        return chunk


def main() -> None:
    from darkmatter.statistical.coding_structure import codon_position_bias

    acc_map = panel_accession_map()
    print(f"{len(acc_map)} panel accessions to look for", flush=True)

    session = requests.Session()
    resp = session.get(ARCHIVE_URL, stream=True, timeout=(15, 60))
    resp.raise_for_status()
    counter = _CountingRaw(resp.raw)

    rows = []
    genomes_done = 0
    t0 = time.time()
    members_seen = 0

    with tarfile.open(fileobj=counter, mode="r|gz") as tar:
        for member in tar:
            members_seen += 1
            elapsed = time.time() - t0
            if genomes_done >= MAX_GENOMES or counter.total >= MAX_BYTES or elapsed >= MAX_SECONDS:
                print(f"stopping: genomes={genomes_done} bytes={counter.total/1e9:.1f}GB "
                      f"elapsed={elapsed:.0f}s (budget reached)", flush=True)
                break
            if not member.isfile() or not member.name.endswith("_genomic.fna.gz"):
                continue
            fname = member.name.rsplit("/", 1)[-1]
            raw_acc = fname[: -len("_genomic.fna.gz")]
            # archive uses e.g. GCA_000147015.1; strip any trailing assembly
            # suffix mismatch isn't expected here, direct lookup:
            base_acc = raw_acc
            if base_acc not in acc_map:
                continue
            panel_acc = acc_map[base_acc]
            coords = gene_coords(panel_acc)
            if not coords:
                continue

            import gzip
            handle = tar.extractfile(member)
            if handle is None:
                continue
            genome_fasta = gzip.decompress(handle.read()).decode("utf-8", errors="ignore")
            contigs = parse_fasta(genome_fasta)

            n_seg = 0
            for contig, genes in coords.items():
                seq = contigs.get(contig)
                if not seq:
                    continue
                for seg in intergenic_segments(seq, genes):
                    cpbb = codon_position_bias(seg)
                    if not np.isnan(cpbb):
                        rows.append({"genome_accession": panel_acc, "contig": contig,
                                    "length": len(seg), "cpbb": cpbb})
                        n_seg += 1
            genomes_done += 1
            print(f"  [{genomes_done}/{MAX_GENOMES}] {panel_acc}: {n_seg} intergenic segments "
                  f"({counter.total/1e9:.2f}GB read, {elapsed:.0f}s)", flush=True)

    if not rows:
        print("no intergenic segments collected -- nothing to write", flush=True)
        return

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}: {len(df)} segments from {df['genome_accession'].nunique()} genomes", flush=True)

    # compare against the real dark-gene CPBB already on disk (Layer 5, P4-D5)
    real = pd.read_csv(PROC / "coding_structure.csv")["cpbb"].dropna().to_numpy()
    intergenic = df["cpbb"].to_numpy()

    def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
        y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
        s = np.r_[pos, neg]
        order = np.argsort(s, kind="mergesort")
        ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s) + 1)
        p = y == 1
        return float((ranks[p].sum() - p.sum() * (p.sum() + 1) / 2) / (p.sum() * (~p).sum()))

    result = {
        "n_genomes_scanned": genomes_done, "n_segments": int(len(df)),
        "bytes_read_gb": counter.total / 1e9,
        "intergenic_cpbb_median": float(np.median(intergenic)),
        "real_cpbb_median": float(np.median(real)),
        "coding_vs_true_intergenic_auroc": auroc(real, intergenic),
        "scope_note": ("bounded prefix of the 65GB genome archive, not the full 502-genome "
                       "panel -- see script docstring; full-panel coverage is optional future work"),
    }
    RESULTS_OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nreal CPBB median {result['real_cpbb_median']:.4f} vs "
          f"true-intergenic median {result['intergenic_cpbb_median']:.4f}")
    print(f"coding-vs-true-intergenic AUROC = {result['coding_vs_true_intergenic_auroc']:.4f}")
    print(f"wrote {RESULTS_OUT}")


if __name__ == "__main__":
    main()
