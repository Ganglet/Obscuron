"""Layer 4 -- multi-signal convergence (blueprint SS6, the Phase-4 extension).

Combines independent novelty axes and assigns high confidence only where they
agree (the SS4 "independent lines converge" principle). Here the two axes are
the Layer-1 ESM-2 EVT embedding novelty and a genomic-context novelty computed
from contig gene-order -- deliberately independent (one is sequence-embedding
distance, the other is genomic neighbourhood), so their agreement is genuine
multi-evidence rather than one signal restated.
"""

from darkmatter.convergence.multisignal import (
    aa_composition,
    convergence_stats,
    genomic_context_novelty,
    high_confidence_convergent,
    high_confidence_convergent_multi,
    mahalanobis_novelty,
)

__all__ = [
    "genomic_context_novelty",
    "convergence_stats",
    "high_confidence_convergent",
    "high_confidence_convergent_multi",
    "aa_composition",
    "mahalanobis_novelty",
]
