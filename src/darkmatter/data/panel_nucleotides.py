"""Extract genome-panel nucleotide CDS FASTA directly from GTDB's public
mirror -- the modality Genos-m actually needs.

Everything fetched so far (S3 and locally) is `gtdb_proteins_aa_reps` --
amino-acid protein sequences, correct for ESM-2 and the Pfam/HMMER
differencing pipeline, but the wrong input for Genos-m, which does
single-nucleotide tokenization (see embeddings/genos_m.py). Protein can't
be reversed back into its source DNA (the genetic code is degenerate --
many codons map to the same amino acid), so this has to come from a
fresh fetch, not a transform of what's already extracted.

Confirmed against the live archive: `gtdb_proteins_nt_reps_r{N}.tar.gz`
exists at the expected GTDB URL, supports Range requests, and has the
identical member layout as the protein archive --
`protein_fna_reps/{domain}/{RS_|GB_}{accession}_protein.fna` -- same
representative genomes, same gene calls, just nucleotide instead of
translated. That means the panel's existing accession list applies
unchanged, and it's the correct like-for-like input for a fair
Genos-m-vs-ESM-2 comparison on the same genes (blueprint D3's own framing
of the asymmetry between the two models).

Streamed straight from the public mirror, not through S3 -- consistent
with Track 1's P1-D4 standard (stream what's needed, don't warehouse),
and it means no AWS credentials are required for this step at all.
"""

from __future__ import annotations

import math
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

from darkmatter.data.panel_proteins import extract_panel_proteins_stream

BASE_URL = "https://data.gtdb.ecogenomic.org/releases"
NT_REPS_FILE = "genomic_files_reps/gtdb_proteins_nt_reps_r{n}.tar.gz"


class _ChunkedHTTPReader:
    """Same bounded-chunk, low-concurrency pattern as _ChunkedS3Reader in
    panel_proteins.py, over plain HTTP Range requests instead of S3
    GetObject -- this connection has repeatedly proven fragile on long
    open streams and high concurrency, regardless of which service is on
    the other end."""

    def __init__(self, url: str, chunk_bytes: int = 16 * 1024 * 1024, max_retries: int = 10, workers: int = 2):
        self._url = url
        self._chunk_bytes = chunk_bytes
        self._max_retries = max_retries
        self._session = requests.Session()

        head = self._session.head(url, timeout=(15, 60))
        head.raise_for_status()
        self._total_bytes = int(head.headers["Content-Length"])
        self._num_chunks = math.ceil(self._total_bytes / chunk_bytes) if self._total_bytes else 0

        self._executor = ThreadPoolExecutor(max_workers=workers)
        self._window = workers * 2
        self._next_to_dispatch = 0
        self._pending: deque = deque()
        self._buf = b""
        self._buf_pos = 0
        self._fill_window()

    def _fetch_chunk(self, idx: int) -> bytes:
        start = idx * self._chunk_bytes
        end = min(start + self._chunk_bytes, self._total_bytes) - 1
        for attempt in range(1, self._max_retries + 1):
            try:
                resp = self._session.get(self._url, headers={"Range": f"bytes={start}-{end}"}, timeout=(15, 60))
                resp.raise_for_status()
                return resp.content
            except Exception as e:
                if attempt == self._max_retries:
                    raise
                wait = min(5 * attempt, 30)
                print(f"  http chunk {idx}/{self._num_chunks} failed ({e}), retrying in {wait}s...", flush=True)
                time.sleep(wait)
        raise RuntimeError("unreachable")

    def _fill_window(self) -> None:
        while len(self._pending) < self._window and self._next_to_dispatch < self._num_chunks:
            idx = self._next_to_dispatch
            self._pending.append(self._executor.submit(self._fetch_chunk, idx))
            self._next_to_dispatch += 1

    def _fetch_next_chunk(self) -> bytes:
        if not self._pending:
            return b""
        data = self._pending.popleft().result()
        self._fill_window()
        return data

    def read(self, amt: int | None = None) -> bytes:
        if amt is None:
            parts = [self._buf[self._buf_pos :]]
            self._buf, self._buf_pos = b"", 0
            while self._pending:
                parts.append(self._fetch_next_chunk())
            return b"".join(parts)

        while len(self._buf) - self._buf_pos < amt and self._pending:
            self._buf = self._buf[self._buf_pos :] + self._fetch_next_chunk()
            self._buf_pos = 0

        data = self._buf[self._buf_pos : self._buf_pos + amt]
        self._buf_pos += len(data)
        return data


def extract_panel_nucleotides_from_gtdb(release: str, accessions: set[str], out_dir: Path) -> list[Path]:
    n = release.lstrip("Rr")
    url = f"{BASE_URL}/release{n}/{n}.0/{NT_REPS_FILE.format(n=n)}"
    stream = _ChunkedHTTPReader(url)
    return extract_panel_proteins_stream(stream, accessions, out_dir, suffix="_protein.fna")
