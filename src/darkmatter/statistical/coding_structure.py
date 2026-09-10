"""Coding-structure statistics for Layer 5.

Two model-free signals of authentic coding structure:
  - codon-position base bias (CPBB): real protein-coding DNA has distinct base
    composition at the three codon positions (e.g. GC3 skew, stop-avoidance in
    frame). A mononucleotide shuffle destroys this; a Markov-1 null preserves
    dinucleotide frequencies but not the reading frame, so it also collapses it.
    CPBB is the frame's chi-square-like divergence from the overall composition.
  - k-mer (n-gram) entropy: Shannon entropy of the amino-acid k-mer distribution;
    real proteins are less uniform than a composition-matched shuffle.

The null controls make this a *discriminator*: real dark genes should score far
above their own shuffled/Markov nulls, evidence they are genuine ORFs and not
spurious calls (the blueprint SS10 non-coding-artifact risk).
"""

from __future__ import annotations

import math
from collections import Counter

import numpy as np

_BASES = "ACGT"


def codon_position_bias(nt: str) -> float:
    """Chi-square-like divergence of the three codon-position base distributions
    from the sequence's overall base distribution. ~0 for frame-less noise, high
    for real coding DNA. Uses only A/C/G/T; length trimmed to a multiple of 3."""
    seq = [c for c in nt.upper() if c in _BASES]
    n = (len(seq) // 3) * 3
    if n < 30:  # too short to estimate three positional distributions
        return float("nan")
    seq = seq[:n]
    overall = Counter(seq)
    tot = sum(overall.values())
    p = {b: overall.get(b, 0) / tot for b in _BASES}
    stat = 0.0
    for j in range(3):
        col = seq[j::3]
        m = len(col)
        cj = Counter(col)
        for b in _BASES:
            if p[b] > 0:
                pbj = cj.get(b, 0) / m
                stat += (pbj - p[b]) ** 2 / p[b]
    return stat / 3.0


def kmer_entropy(seq: str, k: int) -> float:
    """Shannon entropy (bits) of the k-mer distribution of a sequence."""
    if len(seq) < k:
        return float("nan")
    counts = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    tot = sum(counts.values())
    return -sum((c / tot) * math.log2(c / tot) for c in counts.values())


def shuffle_seq(seq: str, rng: np.random.Generator) -> str:
    """Random permutation of the sequence (preserves composition, destroys order)."""
    arr = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    rng.shuffle(arr)
    return arr.tobytes().decode()


def markov1_null(seq: str, rng: np.random.Generator) -> str:
    """Order-1 Markov surrogate: preserves single- and di-symbol frequencies but
    not longer-range / reading-frame structure. A stronger null than a plain shuffle."""
    if len(seq) < 2:
        return seq
    trans: dict[str, list[str]] = {}
    for a, b in zip(seq, seq[1:]):
        trans.setdefault(a, []).append(b)
    out = [seq[0]]
    for _ in range(len(seq) - 1):
        nxt = trans.get(out[-1])
        out.append(nxt[rng.integers(len(nxt))] if nxt else seq[rng.integers(len(seq))])
    return "".join(out)


def coding_structure_features(aa: str, nt: str) -> dict[str, float]:
    """The per-gene Layer-5 feature bundle."""
    return {
        "cpbb": codon_position_bias(nt),
        "aa_entropy_k2": kmer_entropy(aa, 2),
        "aa_entropy_k3": kmer_entropy(aa, 3),
    }
