"""Layer 5 -- statistical discrimination of coding structure from noise (blueprint SS6).

Distinguishes genuinely coding dark-matter genes from spurious / non-coding
artifacts using n-gram entropy and positional (codon) structure statistics, with
shuffled and Markov nulls as the control group (the SS4 computational-linguistics
parallel: statistical regularity within a sequence is evidence of authentic coding
structure). Model-free, no GPU -- pure sequence statistics.
"""

from darkmatter.statistical.coding_structure import (
    codon_position_bias,
    coding_structure_features,
    kmer_entropy,
    markov1_null,
    shuffle_seq,
)

__all__ = [
    "codon_position_bias",
    "kmer_entropy",
    "shuffle_seq",
    "markov1_null",
    "coding_structure_features",
]
