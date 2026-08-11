"""Score panel proteins against a Pfam HMM library at GA (gathering)
threshold — the "dark-at-T0" / "characterised-at-T0" half of P1-D3: a
protein with no Pfam-35 hit at GA is dark at T0; one with a hit is
characterised-at-T0 (the known-boundary reference set).

Runs in the hmmsearch direction (each Pfam HMM against the whole protein
set) rather than hmmscan (each protein against the whole HMM set) — same
GA-threshold semantics, but the natural fast direction in pyhmmer for
"every family against many proteins" rather than a pressed/indexed HMM
database.
"""

from __future__ import annotations

from pathlib import Path

import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile


def load_protein_sequences(faa_path: Path, alphabet: Alphabet) -> DigitalSequenceBlock:
    with SequenceFile(faa_path, digital=True, alphabet=alphabet) as sf:
        return DigitalSequenceBlock(alphabet, list(sf))


def load_protein_sequences_multi(faa_paths: list[Path], alphabet: Alphabet) -> DigitalSequenceBlock:
    """Combine many genomes' proteins into one block so a Pfam library only
    needs to be scanned once across the whole panel, not once per genome.
    Protein IDs are contig-accession-based and unique per genome's own
    assembly, so no collision risk combining across genomes."""
    seqs = []
    for faa_path in faa_paths:
        with SequenceFile(faa_path, digital=True, alphabet=alphabet) as sf:
            seqs.extend(sf)
    return DigitalSequenceBlock(alphabet, seqs)


def scan_against_pfam(
    sequences: DigitalSequenceBlock, pfam_hmm_path: Path, cpus: int = 0
) -> dict[str, list[str]]:
    """Return {sequence_id: [pfam_accessions_hit_at_GA, ...]}, empty list = dark."""
    hits_by_seq: dict[str, list[str]] = {seq.name: [] for seq in sequences}

    with pyhmmer.plan7.HMMFile(pfam_hmm_path) as hmm_file:
        for top_hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="gathering", cpus=cpus):
            pfam_acc = top_hits.query.accession
            for hit in top_hits:
                if hit.included:
                    hits_by_seq[hit.name].append(pfam_acc)

    return hits_by_seq


def scan_genome_against_pfam(faa_path: Path, pfam_hmm_path: Path, cpus: int = 0) -> dict[str, list[str]]:
    alphabet = Alphabet.amino()
    sequences = load_protein_sequences(faa_path, alphabet)
    return scan_against_pfam(sequences, pfam_hmm_path, cpus=cpus)
