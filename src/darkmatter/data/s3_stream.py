"""Stream a URL into S3 with true resume — no local disk buffering, and a
dropped connection continues from the last uploaded byte instead of
restarting the whole file.

Needed for the GTDB representative-genome protein FASTA archives
(R232: ~123GB, R207: ~40GB) — see docs/reproducibility.md "GTDB snapshot
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

Resume works because S3 itself is the source of truth for progress, not a
local state file: on start, list_multipart_uploads finds any in-progress
upload for this key, list_parts says which parts already landed, and the
GTDB fetch resumes from the next byte after that (confirmed their server
supports Range requests — Accept-Ranges: bytes, 206 responses).
"""

from __future__ import annotations

import time

import boto3
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

PART_BYTES = 50 * 1024 * 1024  # smaller than S3's 100MB default: fits well inside whatever
# connection-duration limit is killing long transfers on this network, confirmed by trial.


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


def resumable_stream_to_s3(url: str, bucket: str, key: str, profile: str) -> None:
    s3 = boto3.Session(profile_name=profile).client("s3")
    http = _http_session()

    head = http.head(url, timeout=(15, 60))
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

    completed_parts = [{"PartNumber": p["PartNumber"], "ETag": p["ETag"]} for p in parts]
    bytes_done = sum(p["Size"] for p in parts)
    next_part_number = (parts[-1]["PartNumber"] + 1) if parts else 1

    start_time = time.monotonic()
    last_report = start_time

    while bytes_done < total_bytes:
        range_start = bytes_done
        range_end = min(bytes_done + PART_BYTES, total_bytes) - 1

        data = _fetch_range(http, url, range_start, range_end)

        part_resp = s3.upload_part(Bucket=bucket, Key=key, UploadId=upload_id, PartNumber=next_part_number, Body=data)
        completed_parts.append({"PartNumber": next_part_number, "ETag": part_resp["ETag"]})
        bytes_done += len(data)
        next_part_number += 1

        now = time.monotonic()
        if now - last_report > 30:
            last_report = now
            elapsed = now - start_time
            rate_mb_s = (bytes_done / 1e6) / elapsed if elapsed > 0 else 0.0
            pct = 100 * bytes_done / total_bytes
            eta_min = (total_bytes - bytes_done) / 1e6 / rate_mb_s / 60 if rate_mb_s > 0 else float("inf")
            print(
                f"[{key}] {bytes_done / 1e9:.2f}GB / {total_bytes / 1e9:.2f}GB "
                f"({pct:.1f}%) {rate_mb_s:.1f}MB/s ETA {eta_min:.0f}min",
                flush=True,
            )

    completed_parts.sort(key=lambda p: p["PartNumber"])
    s3.complete_multipart_upload(Bucket=bucket, Key=key, UploadId=upload_id, MultipartUpload={"Parts": completed_parts})
    print(f"[{key}] done", flush=True)
