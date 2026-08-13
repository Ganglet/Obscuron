"""Pfam-37 net-new-family proxy for "characterised by T1".

Angshuman's P1-D3 spec calls for InterPro-latest as the T1 signal, not
Pfam-37 alone — full InterProScan (~6.6GB, multi-database Docker
pipeline) is the correct way to get that, but is a much heavier lift than
anything built so far. This is the interim proxy: a protein counts as
characterised-by-T1 if it hits a Pfam-37 family at GA threshold whose
family ID didn't exist in Pfam-35 — the same net-new-family logic behind
P1-D5's proxy estimate (Pfam 35.0: 19,632 families -> 37.0: 21,979,
+2,347 net new). Narrower than true InterPro-latest: misses cases where
a sequence is only newly characterised via a non-Pfam InterPro member
database. Flag this as a proxy when reporting to Track 1, not as the
final number.
"""

from __future__ import annotations

import gzip
from pathlib import Path


def load_family_ids(clans_tsv_path: Path) -> set[str]:
    """Pfam family IDs (unversioned, e.g. 'PF00001') from Pfam-A.clans.tsv.gz."""
    with gzip.open(clans_tsv_path, "rt", encoding="utf-8") as f:
        return {line.split("\t", 1)[0] for line in f if line.strip()}


def new_families_since(old_clans_tsv: Path, new_clans_tsv: Path) -> set[str]:
    return load_family_ids(new_clans_tsv) - load_family_ids(old_clans_tsv)


def _strip_version(pfam_accession: str) -> str:
    """'PF10417.12' -> 'PF10417'."""
    return pfam_accession.split(".", 1)[0]


def label_protein(pfam35_hits: list[str], pfam37_hits: list[str], new_family_ids: set[str]) -> dict:
    dark_at_t0 = len(pfam35_hits) == 0
    characterised_t1_proxy = any(_strip_version(acc) in new_family_ids for acc in pfam37_hits)
    return {
        "dark_at_t0": dark_at_t0,
        "characterised_t1_proxy": characterised_t1_proxy,
        "positive_proxy": dark_at_t0 and characterised_t1_proxy,
    }
