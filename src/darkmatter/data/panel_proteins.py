"""Extract genome-panel protein FASTA from the combined GTDB rep-protein
tarball, without downloading the whole archive (~40GB for R207, ~123GB
for R232) to disk.

The archives streamed to S3 by s3_stream.py are single tar.gz files, one
member per genome:
`protein_faa_reps/{domain}/{RS_|GB_}{accession}_protein.faa` — member
names match the `accession` column in genome_panel.csv exactly (confirmed
by inspecting the archive directly), so no NCBI assembly_summary lookup
is needed to go from "these 502 accessions" to "these 502 protein files".

Reads the archive as a single sequential gzip+tar stream from S3 and
writes out only the members whose accession is in the requested panel.
gzip decompression is inherently sequential (no random seek), so the
full archive is still read start-to-finish once — but it's never
buffered whole, on S3 egress or on local disk, and disk usage stays
proportional to the panel size, not the release size.
"""

from __future__ import annotations

import tarfile
import time
from pathlib import Path
from typing import BinaryIO

import boto3
from botocore.config import Config


class _ResilientS3Reader:
    """File-like wrapper around S3 GetObject that reconnects with a Range
    request on read failure, resuming from the last byte offset.

    A single sequential read through a 40GB+ archive on a slow/flaky
    connection (same one that needed retry logic for uploads in
    s3_stream.py) will hit read timeouts; without this, one stall restarts
    the whole multi-hour scan from byte 0.
    """

    def __init__(self, s3, bucket: str, key: str, max_retries: int = 10):
        self._s3 = s3
        self._bucket = bucket
        self._key = key
        self._max_retries = max_retries
        self._offset = 0
        self._body = self._open()

    def _open(self):
        kwargs = {"Bucket": self._bucket, "Key": self._key}
        if self._offset:
            kwargs["Range"] = f"bytes={self._offset}-"
        return self._s3.get_object(**kwargs)["Body"]

    def read(self, amt: int | None = None) -> bytes:
        for attempt in range(1, self._max_retries + 1):
            try:
                chunk = self._body.read(amt)
                self._offset += len(chunk)
                return chunk
            except Exception as e:
                if attempt == self._max_retries:
                    raise
                wait = min(5 * attempt, 30)
                print(f"  s3 read failed at byte {self._offset} ({e}), reconnecting in {wait}s...", flush=True)
                time.sleep(wait)
                self._body = self._open()
        raise RuntimeError("unreachable")


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
    stream = _ResilientS3Reader(s3, bucket, key)
    return extract_panel_proteins_stream(stream, accessions, out_dir)
