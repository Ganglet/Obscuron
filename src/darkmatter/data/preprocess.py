"""Parse raw GTDB taxonomy files into a structured table, plus basic counts.

Week 1 scope only: turn `bac120_taxonomy_r{N}.tsv.gz` / `ar53_taxonomy_r{N}.tsv.gz`
(accession + one semicolon-joined lineage string per genome) into a table with
one column per rank, and report how many genomes/taxa are in each snapshot.
Snapshot differencing — matching genomes across the historical/current pair to
find what got characterized in between — is Week 2, once Track 1 fixes the
operational definition of "dark matter at T0".
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

RANK_PREFIXES = {
    "d__": "domain",
    "p__": "phylum",
    "c__": "class_",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}


def _parse_lineage(lineage: str) -> dict[str, str]:
    ranks = {name: "" for name in RANK_PREFIXES.values()}
    for token in lineage.split(";"):
        prefix, _, value = token.partition("__")
        rank = RANK_PREFIXES.get(prefix + "__")
        if rank:
            ranks[rank] = value
    return ranks


def _load_taxonomy_file(path: Path, source: str) -> pd.DataFrame:
    rows = []
    with pd.io.common.get_handle(path, "r", compression="gzip") as handle:
        for line in handle.handle:
            accession, _, lineage = line.strip().partition("\t")
            rows.append({"accession": accession, "source": source, **_parse_lineage(lineage)})
    return pd.DataFrame(rows)


def load_taxonomy_release(release: str, data_root: Path) -> pd.DataFrame:
    """Combined bac120 + ar53 taxonomy table for one GTDB release."""
    n = release.lstrip("Rr")
    release_dir = data_root / f"gtdb_{release}"
    bac = _load_taxonomy_file(release_dir / f"bac120_taxonomy_r{n}.tsv.gz", source="bac120")
    ar = _load_taxonomy_file(release_dir / f"ar53_taxonomy_r{n}.tsv.gz", source="ar53")
    return pd.concat([bac, ar], ignore_index=True)


def summarize(df: pd.DataFrame) -> dict:
    return {
        "total_genomes": len(df),
        "bacteria": int((df["source"] == "bac120").sum()),
        "archaea": int((df["source"] == "ar53").sum()),
        "unique_phyla": int(df["phylum"].nunique()),
        "unique_genera": int(df["genus"].nunique()),
        "unique_species": int(df["species"].nunique()),
        "unclassified_species": int((df["species"] == "").sum()),
    }
