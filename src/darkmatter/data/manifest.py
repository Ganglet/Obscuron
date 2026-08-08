"""Dataset manifest + summary statistics writer.

Track 2 deliverable for Phase 1 (blueprint §8): "Produce the dataset
manifest and preliminary summary statistics." Records what was fetched,
when, and from where, so any downstream result can be traced back to an
exact snapshot — the "documentation and reproducibility standards" Track 1
specifies content for (docs/reproducibility.md) and this writes to.
"""

from __future__ import annotations

import json
import platform
from datetime import datetime, timezone
from pathlib import Path


def write_manifest(out_path: Path, source: str, release: str, files: list[Path], extra: dict | None = None) -> Path:
    entry = {
        "source": source,
        "release": release,
        "fetched_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": [
            {"path": str(f), "size_bytes": f.stat().st_size if f.exists() else None} for f in files
        ],
        "host": platform.node(),
        "extra": extra or {},
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = []
    if out_path.exists():
        manifest = json.loads(out_path.read_text(encoding="utf-8"))
    manifest.append(entry)
    out_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return out_path
