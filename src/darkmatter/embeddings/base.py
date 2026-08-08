"""Common interface both embedding backends implement.

Keeping Genos-m and ESM-2 behind the same interface is what makes them the
"parallel, interchangeable pipelines" the blueprint calls for (§7, §8) —
Phase 1 compares them directly, and either can become the primary source
without touching downstream code.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class Embedder(ABC):
    name: str
    embedding_dim: int

    @abstractmethod
    def embed(self, sequences: list[str], batch_size: int = 8) -> np.ndarray:
        """Return an (len(sequences), embedding_dim) array of mean-pooled embeddings."""
        raise NotImplementedError
