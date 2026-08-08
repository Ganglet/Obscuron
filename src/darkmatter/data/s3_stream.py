"""Stream a URL directly into S3, multipart, without buffering the whole
file to local disk.

Needed for the GTDB representative-genome protein FASTA archives
(R232: ~123GB, R207: ~40GB) — see docs/reproducibility.md "GTDB snapshot
fetch". Neither this machine's C: (near-full) nor E: (85GB free) can hold
either file whole, let alone both, so download and upload happen
concurrently: bytes are read from the HTTP response and handed to boto3's
multipart uploader as they arrive.

Caveat: a network failure partway through is not resumed — the multipart
upload is aborted and the whole transfer must be re-run. Acceptable for a
one-off dataset pull; not built for unreliable links.
"""

from __future__ import annotations

import time

import boto3
import requests
from boto3.s3.transfer import TransferConfig
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

MULTIPART_CHUNK_BYTES = 100 * 1024 * 1024  # 100MB parts -> well under S3's 10,000-part cap


class _StreamAsFile:
    """Adapts a streaming requests.Response into the read(size)-only interface boto3 upload_fileobj needs."""

    def __init__(self, resp: requests.Response, total_bytes: int, label: str):
        self._iter = resp.iter_content(chunk_size=1 << 20)
        self._buf = b""
        self._read_bytes = 0
        self._total_bytes = total_bytes
        self._label = label
        self._last_report = time.monotonic()
        self._start = time.monotonic()

    def read(self, size: int = -1) -> bytes:
        while size < 0 or len(self._buf) < size:
            try:
                chunk = next(self._iter)
            except StopIteration:
                break
            self._buf += chunk
        if size < 0:
            out, self._buf = self._buf, b""
        else:
            out, self._buf = self._buf[:size], self._buf[size:]
        self._read_bytes += len(out)
        self._maybe_report()
        return out

    def _maybe_report(self) -> None:
        now = time.monotonic()
        if now - self._last_report < 30:
            return
        self._last_report = now
        elapsed = now - self._start
        rate_mb_s = (self._read_bytes / 1e6) / elapsed if elapsed > 0 else 0.0
        pct = 100 * self._read_bytes / self._total_bytes if self._total_bytes else 0.0
        eta_min = ((self._total_bytes - self._read_bytes) / 1e6 / rate_mb_s / 60) if rate_mb_s > 0 else float("inf")
        print(
            f"[{self._label}] {self._read_bytes / 1e9:.2f}GB / {self._total_bytes / 1e9:.2f}GB "
            f"({pct:.1f}%) {rate_mb_s:.1f}MB/s ETA {eta_min:.0f}min",
            flush=True,
        )


def stream_to_s3(url: str, bucket: str, key: str, profile: str) -> None:
    session = boto3.Session(profile_name=profile)
    s3 = session.client("s3")

    retry = Retry(total=5, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504])
    http = requests.Session()
    http.mount("https://", HTTPAdapter(max_retries=retry))

    with http.get(url, stream=True, timeout=(10, 120)) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        print(f"Starting {url} -> s3://{bucket}/{key} ({total / 1e9:.2f}GB)", flush=True)
        fileobj = _StreamAsFile(resp, total, key)
        config = TransferConfig(
            multipart_chunksize=MULTIPART_CHUNK_BYTES,
            multipart_threshold=MULTIPART_CHUNK_BYTES,
        )
        s3.upload_fileobj(fileobj, bucket, key, Config=config)
    print(f"Done: s3://{bucket}/{key}", flush=True)
