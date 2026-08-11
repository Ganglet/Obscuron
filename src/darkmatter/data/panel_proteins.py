"""Extract genome-panel protein FASTA from the combined GTDB rep-protein
tarball, without downloading the whole archive (~40GB for R207, ~123GB
for R232) to disk.

The archives streamed to S3 by s3_stream.py are single tar.gz files, one
member per genome:
`protein_faa_reps/{domain}/{RS_|GB_}{accession}_protein.faa` — member
names match the `accession` column in genome_panel.csv exactly (confirmed
by inspecting the archive directly), so no NCBI assembly_summary lookup
is needed to go from "these 502 accessions" to "these 502 protein files".

Reads the archive as a sequential gzip+tar stream, fetched from S3 in
bounded byte-range chunks (not one continuously-open GetObject — that
timed out constantly on this connection), and writes out only the
members whose accession is in the requested panel. gzip decompression is
inherently sequential (no random seek), so the full archive is still
read start-to-finish once — but it's never buffered whole, on S3 egress
or on local disk, and disk usage stays proportional to the panel size,
not the release size.
"""

from __future__ import annotations

import math
import tarfile
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import BinaryIO

import boto3
from botocore.config import Config


class _ChunkedS3Reader:
    """File-like wrapper around S3 GetObject that fetches bounded byte-range
    chunks with a sliding window of concurrent requests, instead of keeping
    one long streaming read open or fetching chunks strictly one at a time.

    A continuously-open StreamingBody read timed out every ~10-100KB on
    this connection, while a plain bounded 5MB range GET finished in under
    10s — the same "one open-ended connection is fragile, bounded ranges
    survive a blip" lesson s3_stream.py already learned the hard way for
    uploads. That saga also found this connection's real per-connection
    cap is only ~1.5-2.7MB/s regardless of chunking, and that concurrent
    workers are what actually raised aggregate throughput — so chunks are
    dispatched several at a time (in order, retried individually on
    failure) rather than one fetch-then-wait per chunk.
    """

    def __init__(
        self,
        s3,
        bucket: str,
        key: str,
        chunk_bytes: int = 16 * 1024 * 1024,
        max_retries: int = 10,
        workers: int = 6,
    ):
        self._s3 = s3
        self._bucket = bucket
        self._key = key
        self._chunk_bytes = chunk_bytes
        self._max_retries = max_retries
        self._total_bytes = s3.head_object(Bucket=bucket, Key=key)["ContentLength"]
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
                resp = self._s3.get_object(Bucket=self._bucket, Key=self._key, Range=f"bytes={start}-{end}")
                return resp["Body"].read()
            except Exception as e:
                if attempt == self._max_retries:
                    raise
                wait = min(5 * attempt, 30)
                print(f"  s3 chunk {idx}/{self._num_chunks} failed ({e}), retrying in {wait}s...", flush=True)
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


def _member_accession(member_name: str) -> str | None:
    fname = member_name.rsplit("/", 1)[-1]
    if not fname.endswith("_protein.faa"):
        return None
    return fname[: -len("_protein.faa")]


def extract_panel_proteins_stream(stream: BinaryIO, accessions: set[str], out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    remaining = set(accessions)
    written: list[Path] = []
    members_seen = 0

    with tarfile.open(fileobj=stream, mode="r|gz") as tar:
        for member in tar:
            members_seen += 1
            if not member.isfile():
                continue
            accession = _member_accession(member.name)
            if accession is None or accession not in remaining:
                continue

            handle = tar.extractfile(member)
            if handle is None:
                continue
            dest = out_dir / f"{accession}_protein.faa"
            dest.write_bytes(handle.read())
            written.append(dest)
            remaining.discard(accession)

            if members_seen % 5000 == 0:
                print(
                    f"  scanned {members_seen} members, found {len(written)}/{len(accessions)} panel genomes so far",
                    flush=True,
                )
            if not remaining:
                break

    if remaining:
        preview = sorted(remaining)[:10]
        suffix = "..." if len(remaining) > 10 else ""
        print(f"  WARNING: {len(remaining)} panel accessions not found in archive: {preview}{suffix}", flush=True)
    return written


def extract_panel_proteins_from_s3(bucket: str, key: str, profile: str, accessions: set[str], out_dir: Path) -> list[Path]:
    config = Config(read_timeout=60, connect_timeout=15, retries={"max_attempts": 3})
    s3 = boto3.Session(profile_name=profile).client("s3", config=config)
    stream = _ChunkedS3Reader(s3, bucket, key)
    return extract_panel_proteins_stream(stream, accessions, out_dir)
