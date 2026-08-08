"""Generic Pfam-A release fetcher — parameterized by release version.

URL layout confirmed against https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/ :
    Pfam-A.clans.tsv.gz   (family -> clan/description list, small, few MB)
    Pfam-A.hmm.gz         (full HMM library, large, hundreds of MB+)

`metadata_only=True` fetches just Pfam-A.clans.tsv.gz — the family list
needed to diff "which families exist in release X vs release Y", without
pulling the full HMM library (only needed later for actually running
hmmscan against historical vs current families).
"""

from __future__ import annotations

from pathlib import Path

from darkmatter.data.download import download

BASE_URL = "http://ftp.ebi.ac.uk/pub/databases/Pfam/releases"

METADATA_FILE = "Pfam-A.clans.tsv.gz"
HMM_FILE = "Pfam-A.hmm.gz"


def fetch_pfam_release(release: str, out_dir: Path, metadata_only: bool = True) -> list[Path]:
    base = f"{BASE_URL}/Pfam{release}"
    filenames = [METADATA_FILE] if metadata_only else [METADATA_FILE, HMM_FILE]

    fetched = []
    for fname in filenames:
        dest = out_dir / f"pfam_{release}" / fname
        if dest.exists():
            fetched.append(dest)
            continue
        download(f"{base}/{fname}", dest)
        fetched.append(dest)
    return fetched
