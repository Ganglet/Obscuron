"""Phylum-stratified genome panel sampler for the Phase 1 retrospective
benchmark (P1-D6, docs/Track1_phase1_benchmark_scope.md section 6).

The full GTDB rep-protein set (~10^8 proteins for R207) can't be embedded
locally and isn't warehoused (P1-D4), so the benchmark draws from a
principled subset: GTDB species-representative genomes, sampled per
phylum with a floor (so rare candidate phyla like CPR/DPANN — where dark
matter concentrates — always appear) and a cap (so Proteobacteria/
Firmicutes don't dominate the panel).

"Species-representative" is a metadata-only fact: GTDB's taxonomy file
(bac120_taxonomy_r{N}.tsv.gz) lists every genome in the release, not just
species reps. The `gtdb_representative` flag that distinguishes them only
exists in the metadata file (bac120_metadata_r{N}), so panel sampling
needs both files, not taxonomy.csv alone.
"""

from __future__ import annotations

import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

_METADATA_BASENAMES = ["bac120_metadata_r{n}", "ar53_metadata_r{n}"]
_METADATA_EXTS = [".tsv.gz", ".tar.gz"]
_METADATA_COLUMNS = ["accession", "gtdb_representative"]


def _read_metadata_columns(path: Path, columns: list[str]) -> pd.DataFrame:
    if path.name.endswith(".tar.gz"):
        with tarfile.open(path, "r:gz") as tar:
            member = next(m for m in tar.getmembers() if m.name.endswith(".tsv"))
            handle = tar.extractfile(member)
            return pd.read_csv(handle, sep="\t", usecols=columns, low_memory=False)
    return pd.read_csv(path, sep="\t", compression="gzip", usecols=columns, low_memory=False)


def load_representative_accessions(release: str, data_root: Path) -> set[str]:
    """Accessions flagged `gtdb_representative == 't'` across bac120 + ar53 metadata."""
    n = release.lstrip("Rr")
    release_dir = data_root / f"gtdb_{release}"

    accessions: set[str] = set()
    for basename_template in _METADATA_BASENAMES:
        basename = basename_template.format(n=n)
        path = next(
            (release_dir / f"{basename}{ext}" for ext in _METADATA_EXTS if (release_dir / f"{basename}{ext}").exists()),
            None,
        )
        if path is None:
            raise FileNotFoundError(f"no metadata file for {basename} under {release_dir} (tried {_METADATA_EXTS})")
        df = _read_metadata_columns(path, _METADATA_COLUMNS)
        accessions.update(df.loc[df["gtdb_representative"] == "t", "accession"])
    return accessions


def allocate_panel_sizes(
    phylum_counts: dict[str, int], target_total: int = 500, floor: int = 2, cap: int = 30, tol: int = 5, max_iter: int = 60
) -> tuple[dict[str, int], float]:
    """Proportional-with-floor-cap allocation (P1-D6): binary-search the
    sampling fraction `f` so per-phylum sizes clip(round(f*n), floor, cap)
    sum to ~target_total. Monotonic in f, so binary search converges."""

    def sizes_for(f: float) -> dict[str, int]:
        return {phylum: min(n, max(floor, min(cap, round(f * n)))) for phylum, n in phylum_counts.items()}

    lo, hi, f = 0.0, 1.0, 1.0
    for _ in range(max_iter):
        f = (lo + hi) / 2
        total = sum(sizes_for(f).values())
        if abs(total - target_total) <= tol:
            break
        if total < target_total:
            lo = f
        else:
            hi = f
    return sizes_for(f), f


def sample_panel(
    taxonomy: pd.DataFrame,
    representative_accessions: set[str],
    target_total: int = 500,
    floor: int = 2,
    cap: int = 30,
    seed: int = 42,
) -> tuple[pd.DataFrame, dict]:
    reps = taxonomy[taxonomy["accession"].isin(representative_accessions)].drop_duplicates("accession")
    phylum_counts = reps.groupby("phylum")["accession"].nunique().to_dict()
    sizes, f = allocate_panel_sizes(phylum_counts, target_total=target_total, floor=floor, cap=cap)

    rng = np.random.default_rng(seed)
    chosen = [
        reps[reps["phylum"] == phylum].sample(n=n, random_state=rng)
        for phylum, n in sizes.items()
        if n > 0
    ]
    panel = pd.concat(chosen, ignore_index=True)

    info = {
        "target_total": target_total,
        "floor": floor,
        "cap": cap,
        "seed": seed,
        "allocation_f": f,
        "panel_size": len(panel),
        "phyla_available": len(phylum_counts),
        "phyla_in_panel": panel["phylum"].nunique(),
        "sizes_by_phylum": sizes,
    }
    return panel, info
