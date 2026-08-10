"""Stream a URL into S3 with true resume — no local disk buffering, and a
dropped connection continues from the last uploaded byte instead of
restarting the whole file.

Needed for the GTDB representative-genome protein FASTA archives
(R232: ~132GB, R207: ~43GB) — see docs/reproducibility.md "GTDB snapshot
fetch". Neither this machine's C: (near-full) nor E: (85GB free) can hold
either file whole, so download and upload happen concurrently: each part
is fetched as its own bounded HTTP Range request and handed straight to
an S3 multipart part — never the whole file in memory, never on disk.

Each part is its own request, not one long-lived connection for the whole
remaining file. Found necessary in practice: the first live run died 8.5%
into a 43GB file on a plain read timeout; after switching to a resumable
design, a live run still failed repeatedly at an oddly consistent
~100-108MB into every attempt — a strong sign something on the network
path (router/ISP connection tracking, most likely) kills long-lived
connections after a roughly fixed duration, not a fixed byte count. A
single open-ended `Range: bytes=X-` GET streaming tens of GB is doomed to
hit that wall every time no matter how the upload side is chunked.
Requesting `Range: bytes=start-end` per part instead gets a fresh
connection every ~50MB, so a mid-transfer network blip only costs one
part's retry, not the entire rest of the file.

Parts fetch+upload concurrently (small thread pool), not one at a time.
Found necessary after a real wifi-plan upgrade (200Mbps) didn't move
sustained throughput at all — the old sequential "fetch, then upload,
then fetch the next one" design left the pipe idle half the time no
matter how fast it was, plus every part pays a full fresh TLS handshake
(compounded by local SSL-inspection software), so overhead dominated
over actual transfer time. Concurrency overlaps that idle/handshake time
across parts instead of paying it serially.

Resume works because S3 itself is the source of truth for progress, not a
local state file: on start, list_multipart_uploads finds any in-progress
upload for this key, list_parts says which parts already landed, and the
GTDB fetch resumes from the next byte after that (confirmed their server
supports Range requests — Accept-Ranges: bytes, 206 responses). This
holds regardless of how the remaining parts get fetched, so switching to
concurrent workers never risks already-uploaded parts.
"""

from __future__ import annotations

import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import boto3
import requests
from botocore.config import Config
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Recent botocore versions default to attaching a CRC32 checksum requirement
# to every S3 multipart upload (request_checksum_calculation="when_supported").
# upload_part() here doesn't compute/send one, so complete_multipart_upload
# fails at the very end with "missing checksum for part N" — found the hard
# way after a full 43GB transfer completed and then couldn't be finalized.
# "when_required" opts back out: no checksum unless the operation demands one.
_S3_CONFIG = Config(request_checksum_calculation="when_required")

PART_BYTES = 50 * 1024 * 1024  # smaller than S3's 100MB default: fits well inside whatever
# connection-duration limit is killing long transfers on this network, confirmed by trial.

# Conservative on purpose: this network is sensitive to connection churn
# (periodic SSL-inspection failures, a connection-duration limit that
# forced the per-part design above) — too much concurrency risks making
# both worse, not just faster. Dropped from 4 to 2 after a live run showed
# 3 workers hit the identical failure (same byte offset) at the same
# instant — the local SSL-inspection layer appears to interrupt every open
# connection at once periodically, so more concurrent connections just
# means more of them wedge simultaneously, not more throughput.
MAX_WORKERS = 2


def _http_session() -> requests.Session:
    retry = Retry(total=3, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def _find_existing_upload(s3, bucket: str, key: str) -> str | None:
    resp = s3.list_multipart_uploads(Bucket=bucket, Prefix=key)
    for upload in resp.get("Uploads", []):
        if upload["Key"] == key:
            return upload["UploadId"]
    return None


def _existing_parts(s3, bucket: str, key: str, upload_id: str) -> list[dict]:
    parts: list[dict] = []
    marker = None
    while True:
        kwargs = {"Bucket": bucket, "Key": key, "UploadId": upload_id}
        if marker:
            kwargs["PartNumberMarker"] = marker
        resp = s3.list_parts(**kwargs)
        parts.extend(resp.get("Parts", []))
        if not resp.get("IsTruncated"):
            break
        marker = resp.get("NextPartNumberMarker")
    parts.sort(key=lambda p: p["PartNumber"])
    return parts


def _fetch_range(http: requests.Session, url: str, start: int, end: int, max_retries: int = 15) -> bytes:
    """Fetch one bounded byte range, retrying just this range on failure."""
    for attempt in range(1, max_retries + 1):
        try:
            resp = http.get(url, headers={"Range": f"bytes={start}-{end}"}, timeout=(15, 90))
            resp.raise_for_status()
            return resp.content
        except Exception as e:
            if attempt == max_retries:
                raise
            wait = min(5 * attempt, 60)
            print(f"  range {start}-{end} attempt {attempt}/{max_retries} failed ({e}), retrying in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError("unreachable")


def _upload_part(s3, bucket: str, key: str, upload_id: str, part_number: int, data: bytes, max_retries: int = 15) -> dict:
    """upload_part with its own retry — botocore's built-in retries don't
    cover every transient failure (an EndpointConnectionError crashed a live
    run with no retry at all before this existed)."""
    for attempt in range(1, max_retries + 1):
        try:
            return s3.upload_part(Bucket=bucket, Key=key, UploadId=upload_id, PartNumber=part_number, Body=data)
        except Exception as e:
            if attempt == max_retries:
                raise
            wait = min(5 * attempt, 60)
            print(f"  upload_part {part_number} attempt {attempt}/{max_retries} failed ({e}), retrying in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError("unreachable")


def resumable_stream_to_s3(url: str, bucket: str, key: str, profile: str, max_workers: int = MAX_WORKERS) -> None:
    s3 = boto3.Session(profile_name=profile).client("s3", config=_S3_CONFIG)  # boto3 clients are thread-safe

    head = requests.head(url, timeout=(15, 60))
    head.raise_for_status()
    total_bytes = int(head.headers["Content-Length"])

    upload_id = _find_existing_upload(s3, bucket, key)
    if upload_id:
        parts = _existing_parts(s3, bucket, key, upload_id)
        print(f"[{key}] resuming, {len(parts)} parts already uploaded", flush=True)
    else:
        upload_id = s3.create_multipart_upload(Bucket=bucket, Key=key)["UploadId"]
        parts = []
        print(f"[{key}] starting new upload ({total_bytes / 1e9:.2f}GB)", flush=True)

    completed_parts: dict[int, str] = {p["PartNumber"]: p["ETag"] for p in parts}
    bytes_done = sum(p["Size"] for p in parts)
    next_part_number = (parts[-1]["PartNumber"] + 1) if parts else 1

    # Pre-compute every remaining (part_number, start, end) triple up front so
    # part numbers stay deterministic regardless of which worker finishes first.
    ranges = []
    offset = bytes_done
    part_number = next_part_number
    while offset < total_bytes:
        end = min(offset + PART_BYTES, total_bytes) - 1
        ranges.append((part_number, offset, end))
        offset = end + 1
        part_number += 1

    start_time = time.monotonic()
    state_lock = threading.Lock()
    # rate is measured incrementally (since the last report), not cumulatively
    # since this run started — cumulative-since-start divides ALL bytes ever
    # done (including whatever was already uploaded before this resume) by
    # only this run's elapsed time, which massively overstates speed right
    # after every resume and never fully corrects. Track a report_bytes/
    # report_time baseline and reset it each time we print.
    state = {
        "bytes_done": bytes_done,
        "last_report": start_time,
        "report_bytes": bytes_done,
    }

    def process_one(pn: int, start: int, end: int) -> None:
        http = _http_session()  # own session per task, no cross-thread sharing
        data = _fetch_range(http, url, start, end)
        part_resp = _upload_part(s3, bucket, key, upload_id, pn, data)
        with state_lock:
            completed_parts[pn] = part_resp["ETag"]
            state["bytes_done"] += len(data)
            now = time.monotonic()
            interval = now - state["last_report"]
            if interval > 30:
                interval_bytes = state["bytes_done"] - state["report_bytes"]
                rate_mb_s = (interval_bytes / 1e6) / interval if interval > 0 else 0.0
                state["last_report"] = now
                state["report_bytes"] = state["bytes_done"]
                pct = 100 * state["bytes_done"] / total_bytes
                eta_min = (total_bytes - state["bytes_done"]) / 1e6 / rate_mb_s / 60 if rate_mb_s > 0 else float("inf")
                print(
                    f"[{key}] {state['bytes_done'] / 1e9:.2f}GB / {total_bytes / 1e9:.2f}GB "
                    f"({pct:.1f}%) {rate_mb_s:.1f}MB/s ETA {eta_min:.0f}min",
                    flush=True,
                )

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_one, pn, s, e) for pn, s, e in ranges]
        for future in as_completed(futures):
            future.result()  # re-raise if a part exhausted its own retries

    sorted_parts = [{"PartNumber": pn, "ETag": completed_parts[pn]} for pn in sorted(completed_parts)]
    s3.complete_multipart_upload(Bucket=bucket, Key=key, UploadId=upload_id, MultipartUpload={"Parts": sorted_parts})
    print(f"[{key}] done", flush=True)
