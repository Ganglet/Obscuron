"""Slow test: actually loads model weights and runs inference.

Skipped by default (downloads GBs of weights). Run explicitly with:
    RUN_SLOW_TESTS=1 uv run pytest tests/test_embeddings_smoke.py -v
"""

import os

import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("RUN_SLOW_TESTS") != "1",
    reason="downloads model weights; set RUN_SLOW_TESTS=1 to run",
)


def test_esm2_embeds_toy_sequences():
    from darkmatter.embeddings import load_embedder

    embedder = load_embedder("esm2")
    vecs = embedder.embed(["MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNL"])
    assert vecs.shape == (1, embedder.embedding_dim)


def test_genos_m_embeds_toy_sequences():
    from darkmatter.embeddings import load_embedder

    embedder = load_embedder("genos-m")
    vecs = embedder.embed(["ATGCGATCGTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCTGATCGATCGATCGATCG"])
    assert vecs.shape == (1, embedder.embedding_dim)
