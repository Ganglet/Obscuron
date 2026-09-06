"""Layer 3 — immune self/non-self discrimination via calibrated negative selection.

Self = characterised reference; a V-detector repertoire covers the complement
(non-self); the decision threshold is calibrated to a self-tolerance bound.
See docs/Track1_phase3_immune_design.md (P3-D1..D5).
"""

from darkmatter.immune.negative_selection import VDetector

__all__ = ["VDetector"]
