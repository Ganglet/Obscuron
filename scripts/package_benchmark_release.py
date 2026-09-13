#!/usr/bin/env python
"""Package the retrospective dark-matter benchmark (genome panel + protein
labels) as a standalone, reusable release -- independent of the rest of the
Obscuron codebase, so someone who only wants the benchmark (not the scorer,
not the extension layers) can take exactly that.

Ships only derived labels (dark-at-T0 / characterised-by-T1 booleans per
protein, and the panel's genome/taxonomy list), never raw sequences --
consistent with this project's data-provenance standard (P1-D4): sequences
stay streamable from GTDB's public mirror by accession, never re-warehoused.

    uv run python scripts/package_benchmark_release.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "gtdb_R207"
OUT = ROOT / "benchmark_release"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT).strip()
    except Exception:
        return "unknown"


def main() -> None:
    OUT.mkdir(exist_ok=True)

    panel_src = PROC / "genome_panel.csv"
    labels_src = PROC / "panel_protein_labels.csv"
    panel_dst = OUT / "genome_panel.csv"
    labels_dst = OUT / "protein_labels.csv"
    panel_dst.write_bytes(panel_src.read_bytes())
    labels_dst.write_bytes(labels_src.read_bytes())

    panel = pd.read_csv(panel_dst)
    labels = pd.read_csv(labels_dst)
    n_dark = int(labels["dark_at_t0"].sum())
    n_char = int(labels["characterised_t1_proxy"].sum())
    n_pos = int(labels["positive_proxy"].sum())

    manifest = {
        "name": "Obscuron retrospective dark-matter benchmark",
        "version": "1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_commit": _git_commit(),
        "boundary": {
            "T0": {"gtdb_release": "R207", "date": "2022-04", "pfam_release": "35.0", "pfam_date": "2021-11"},
            "T1_proxy": {"gtdb_release": "R232", "date": "2026-04", "pfam_release": "37.0", "pfam_date": "2024-06",
                        "note": "Pfam 37.0 net-new-family proxy stands in for full InterPro-latest "
                                "(P1-D8); narrower than the true signal, so the true positive count "
                                "is a floor, not an exact figure."},
        },
        "sampling": {
            "scheme": "proportional-with-floor-cap by GTDB phylum",
            "target_genomes": 500, "floor": 2, "cap": 30, "seed": 42,
            "actual_panel_size": int(len(panel)), "phyla_covered": int(panel["phylum"].nunique()),
        },
        "labels": {
            "definition": {
                "dark_at_t0": "no Pfam-35.0 domain hit at the family gathering (GA) threshold",
                "characterised_t1_proxy": "hits a Pfam-37.0 family that is new since Pfam-35.0",
                "positive_proxy": "dark_at_t0 AND characterised_t1_proxy",
            },
            "n_proteins": int(len(labels)), "n_dark_at_t0": n_dark,
            "n_characterised_t1_proxy": n_char, "n_positive_proxy": n_pos,
            "dark_fraction": n_dark / len(labels), "positive_fraction_of_all": n_pos / len(labels),
            "go_no_go_floor": "50-100 positives (pre-registered)",
            "margin_over_floor": f"{n_pos / 100:.0f}-{n_pos / 50:.0f}x (against the 50-100 floor range)",
        },
        "files": {
            "genome_panel.csv": {"rows": int(len(panel)), "sha256": _sha256(panel_dst)},
            "protein_labels.csv": {"rows": int(len(labels)), "sha256": _sha256(labels_dst)},
        },
        "not_included": (
            "Raw sequences (protein/nucleotide FASTA) and embeddings are not "
            "shipped: they are large, regenerable, and freely re-downloadable "
            "from GTDB's own public mirror by the accessions in genome_panel.csv. "
            "See README.md in this directory for the exact fetch commands."
        ),
        "license": "MIT (see LICENSE)",
        "citation": "See the parent repository's README.md for the paper citation once published.",
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (OUT / "LICENSE").write_bytes((ROOT / "LICENSE").read_bytes())

    readme = f"""# Obscuron Retrospective Dark-Matter Benchmark

A standalone release of the benchmark behind *Obscuron: Calibrated Novelty
Detection for Microbial Dark Matter via Retrospective Validation*. This
directory is self-contained: it does not require the rest of the Obscuron
repository to be useful on its own.

## What this is

A phylum-stratified panel of {manifest['sampling']['actual_panel_size']} GTDB R207 representative genomes
({manifest['sampling']['phyla_covered']} phyla), with every one of their
{manifest['labels']['n_proteins']:,} proteins labelled against a retrospective
temporal boundary:

- **T0** (historical): GTDB R207 (April 2022) + Pfam 35.0 (November 2021)
- **T1 proxy** (current): GTDB R232 (April 2026) + a Pfam 37.0 net-new-family
  proxy standing in for full InterPro-latest (narrower than the true signal;
  see `manifest.json` for the caveat)

A protein is **dark at T0** if it has no Pfam-35.0 hit at the family
gathering threshold. It is a **positive** if it is dark at T0 *and* matches
a family that is new as of Pfam 37.0, i.e. it went from unannotated to
annotated across the boundary. This is the retrospective ground truth: score
what was dark matter at T0, then check which of it got characterised by T1.

**Counts:** {manifest['labels']['n_dark_at_t0']:,} dark-at-T0 ({manifest['labels']['dark_fraction']*100:.1f}%), {manifest['labels']['n_positive_proxy']:,} positives ({manifest['labels']['positive_fraction_of_all']*100:.2f}% of all proteins) --
{manifest['labels']['margin_over_floor']} the pre-registered go/no-go floor of 50-100 positives.

## Files

| File | Contents |
|---|---|
| `genome_panel.csv` | The {manifest['sampling']['actual_panel_size']} panel genomes: accession, domain, and full GTDB taxonomy (phylum through species). |
| `protein_labels.csv` | One row per protein: `genome_accession`, `protein_id`, `dark_at_t0`, `characterised_t1_proxy`, `positive_proxy`. |
| `manifest.json` | Machine-readable provenance: boundary definition, sampling scheme, label definitions, file checksums. |

## What is *not* included, and why

Raw sequences and embeddings are not shipped. They are large, regenerable,
and freely available from GTDB's own public mirror -- shipping them would
duplicate free data rather than add anything. To fetch the sequences for
this exact panel:

```bash
# protein (amino-acid) sequences for the panel's genomes, streamed directly
# from the public GTDB mirror -- no AWS credentials needed
uv run python scripts/extract_panel_proteins.py --release R207 \\
    --accessions genome_panel.csv

# nucleotide sequences (needed only for a genomic, not protein, embedding model)
uv run python scripts/extract_panel_nucleotides.py --release R207 \\
    --accessions genome_panel.csv
```

(these scripts live in the parent Obscuron repository, not in this release
directory)

## Using the labels directly

```python
import pandas as pd

labels = pd.read_csv("protein_labels.csv")
dark_queries = labels[labels["dark_at_t0"]]
positives = labels[labels["positive_proxy"]]
reference = labels[~labels["dark_at_t0"]]  # characterised-at-T0
```

## Reproducibility

Every checksum, count, and design parameter in this release is recorded in
`manifest.json`, generated by `scripts/package_benchmark_release.py` in the
parent repository at commit `{manifest['source_commit']}`. The full,
numbered decision log behind every choice in this benchmark's construction
(snapshot boundary rationale, sampling scheme, go/no-go gate) is in the
parent repository's `docs/problems_and_decisions.md`.

## License

MIT -- see `LICENSE`.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")

    print(f"wrote {OUT}/ (README.md, manifest.json, LICENSE, genome_panel.csv, protein_labels.csv)")
    print(f"  {manifest['labels']['n_proteins']:,} proteins, {n_dark:,} dark-at-T0, {n_pos:,} positives")


if __name__ == "__main__":
    main()
