"""Shared HTTP download-with-progress helper used by gtdb.py and pfam.py."""

from __future__ import annotations

import time
from pathlib import Path

import requests
from tqdm import tqdm

# A mid-stream connection drop can end resp.iter_content()'s generator
# early without raising -- found the hard way when a 293MB-expected Pfam
# HMM fetch silently landed at 230MB with no exception, and the caller
# had no way to know. Retry the whole download (bounded file, not worth
# resuming byte-by-byte) whenever bytes written don't match Content-Length.


def download(url: str, dest: Path, max_retries: int = 5) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, max_retries + 1):
        written = 0
        with requests.get(url, stream=True, timeout=(10, 120)) as resp:
            resp.raise_for_status()
            total = int(resp.headers.get("content-length", 0))
            with open(dest, "wb") as f, tqdm(total=total, unit="B", unit_scale=True, desc=dest.name) as bar:
                for chunk in resp.iter_content(chunk_size=1 << 16):
                    f.write(chunk)
                    written += len(chunk)
                    bar.update(len(chunk))

        if total == 0 or written == total:
            return

        if attempt == max_retries:
            raise IOError(f"download of {url} truncated after {max_retries} attempts: got {written}/{total} bytes")
        wait = min(5 * attempt, 30)
        print(f"  {dest.name}: got {written}/{total} bytes, retrying in {wait}s...", flush=True)
        time.sleep(wait)
