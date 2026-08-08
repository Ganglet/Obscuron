"""Stream a URL into S3 with true resume — no local disk buffering, and a
dropped connection continues from the last uploaded byte instead of
restarting the whole file.

Needed for the GTDB representative-genome protein FASTA archives
(R232: ~123GB, R207: ~40GB) — see docs/reproducibility.md "GTDB snapshot
fetch". Neither this machine's C: (near-full) nor E: (85GB free) can hold
either file whole, so download and upload happen concurrently: bytes are
read from the HTTP response and handed straight to an S3 multipart part.

Resume works because S3 itself is the source of truth for progress, not a
local state file: on start, list_multipart_uploads finds any in-progress
upload for this key, list_parts says which byte ranges are already
uploaded, and the GTDB fetch resumes from there with a Range request
(confirmed their server supports it — Accept-Ranges: bytes, 206 on a
range GET). A local process restart, laptop sleep, or wifi drop only
costs whatever part was mid-flight when it happened, not the whole file.
Found necessary in practice: the first live run died 8.5% into a 43GB
file on a plain read timeout with no resume, and had to restart from 0.
"""

from __future__ import annotations

import time

import boto3
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

MULTIPART_CHUNK_BYTES = 100 * 1024 * 1024  # 100MB parts -> well under S3's 10,000-part cap


def _http_session() -> requests.Session:
    retry = Retry(total=5, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504])
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


def resumable_stream_to_s3(url: str, bucket: str, key: str, profile: str, max_retries: int = 50) -> None:
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

    for attempt in range(1, max_retries + 1):
        if bytes_done >= total_bytes:
            break
        try:
            with http.get(url, headers={"Range": f"bytes={bytes_done}-"}, stream=True, timeout=(15, 300)) as resp:
                resp.raise_for_status()
                chunks = resp.iter_content(chunk_size=1 << 20)
                buf = bytearray()
                exhausted = False
                while True:
                    while len(buf) < MULTIPART_CHUNK_BYTES and not exhausted:
                        try:
                            buf += next(chunks)
                        except StopIteration:
                            exhausted = True
                    if not buf:
                        break
                    part_bytes = bytes(buf[:MULTIPART_CHUNK_BYTES])
                    del buf[:MULTIPART_CHUNK_BYTES]

                    part_resp = s3.upload_part(
                        Bucket=bucket, Key=key, UploadId=upload_id, PartNumber=next_part_number, Body=part_bytes
                    )
                    completed_parts.append({"PartNumber": next_part_number, "ETag": part_resp["ETag"]})
                    bytes_done += len(part_bytes)
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
            break
        except Exception as e:
            print(f"[{key}] attempt {attempt}/{max_retries} failed at {bytes_done / 1e9:.2f}GB: {e}", flush=True)
            if attempt == max_retries:
                raise
            wait = min(60 * attempt, 600)
            print(f"[{key}] retrying in {wait}s, resuming from {bytes_done / 1e9:.2f}GB...", flush=True)
            time.sleep(wait)

    completed_parts.sort(key=lambda p: p["PartNumber"])
    s3.complete_multipart_upload(Bucket=bucket, Key=key, UploadId=upload_id, MultipartUpload={"Parts": completed_parts})
    print(f"[{key}] done", flush=True)
