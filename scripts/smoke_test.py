#!/usr/bin/env python
"""Loads both embedding backends on this machine and embeds a few toy
sequences, reporting peak VRAM used. This is the real hardware-fit check
for the plan's claim that Genos-m fits in 8-bit on an 8GB-class GPU
(blueprint §9) — not a mock.

Usage:
    uv run python scripts/smoke_test.py
    uv run python scripts/smoke_test.py --skip-genos-m   # ESM-2 only, faster
"""

from __future__ import annotations

import argparse
import time

import torch

from darkmatter.device import detect_device, peak_vram_gb
from darkmatter.embeddings import load_embedder

TOY_DNA = [
    "ATGCGATCGTAGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCTGATCGATCGATCGATCG",
    "GGGGCCCCATATGCGCGCATATATATCGCGCGATATATCGCGCGCGATATATATCGCGCG",
    "TTTTAAAACCCCGGGGATCGATCGATCGATCGATCGTAGCTAGCTAGCTGATCGATCGAT",
]

TOY_PROTEIN = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ",
    "MSEQVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKV",
    "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLLVKHSQ",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-genos-m", action="store_true")
    parser.add_argument("--skip-esm2", action="store_true")
    args = parser.parse_args()

    device, total_vram = detect_device()
    print(f"Detected device: {device}" + (f" ({total_vram:.1f} GB VRAM)" if total_vram else ""))

    if not args.skip_esm2:
        print("\n--- ESM-2 ---")
        t0 = time.time()
        esm2 = load_embedder("esm2")
        print(f"Loaded in {time.time() - t0:.1f}s (dim={esm2.embedding_dim}, quantize_bits={esm2.plan.quantize_bits})")
        vecs = esm2.embed(TOY_PROTEIN)
        print(f"Embedded {len(TOY_PROTEIN)} toy protein sequences -> {vecs.shape}")
        vram = peak_vram_gb()
        if vram is not None:
            print(f"Peak VRAM after ESM-2: {vram:.2f} GB")
        del esm2
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()

    if not args.skip_genos_m:
        print("\n--- Genos-m ---")
        t0 = time.time()
        genos = load_embedder("genos-m")
        print(f"Loaded in {time.time() - t0:.1f}s (dim={genos.embedding_dim}, quantize_bits={genos.plan.quantize_bits})")
        vecs = genos.embed(TOY_DNA)
        print(f"Embedded {len(TOY_DNA)} toy DNA sequences -> {vecs.shape}")
        vram = peak_vram_gb()
        if vram is not None:
            print(f"Peak VRAM after Genos-m: {vram:.2f} GB")

    print("\nSmoke test complete.")


if __name__ == "__main__":
    main()
