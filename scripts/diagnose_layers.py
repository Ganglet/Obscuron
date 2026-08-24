#!/usr/bin/env python
"""Per-layer x centering separation diagnostic, for either arm.

Generalises diagnose_genos_m_layers.py to both models so ESM-2 and Genos-m
are compared on the SAME 100-seq matched set, SAME analysis (raw vs
mean-centered cosine gap, every hidden layer). Reuses the family labels
from genos-m_separation_full.json -> no Pfam rescan. Reads the right
modality per model (ESM-2 <- panel_proteins/*.faa, Genos-m <-
panel_nucleotides/*.fna; same record ids, P1-D9).

    uv run python scripts/diagnose_layers.py --model esm2
    uv run python scripts/diagnose_layers.py --model genos-m   # (already saved)

Answers the P1-D13 layer ablation and whether Genos-m's raw last-layer
collapse (0.0042) was anisotropy/layer choice rather than a real weakness.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

PROC = Path(__file__).resolve().parents[1] / "data" / "processed" / "gtdb_R207"
MODALITY = {  # model -> (sequence dir, glob, needs Genos-m N-run preprocessing)
    "esm2": ("panel_proteins", "*.faa", False),
    "genos-m": ("panel_nucleotides", "*.fna", True),
}


def separation_gap(vectors: np.ndarray, labels: np.ndarray, center: bool) -> float:
    v = vectors - vectors.mean(axis=0, keepdims=True) if center else vectors
    v = v / np.clip(np.linalg.norm(v, axis=1, keepdims=True), 1e-12, None)
    sim = v @ v.T
    same = labels[:, None] == labels[None, :]
    off = ~np.eye(len(labels), dtype=bool)
    return float(sim[same & off].mean() - sim[~same & off].mean())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", choices=["esm2", "genos-m"], required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--fp32", action="store_true",
                    help="cast the model to fp32 -- ESM-2's intermediate hidden states "
                         "overflow fp16 on MPS (NaN layers); needed for a clean layer sweep")
    args = ap.parse_args()
    seq_dir, glob, needs_pre = MODALITY[args.model]

    # 1. exact 100-seq labeled subset from the prior run (no rescan)
    families = json.loads((PROC / "genos-m_separation_full.json").read_text())["families"]
    ordered_ids = [sid for ids in families.values() for sid in ids]
    labels = np.array([fam for fam, ids in families.items() for _ in ids])
    print(f"{len(ordered_ids)} sequences, {len(families)} families (reused)", flush=True)

    # 2. sequences for this modality
    from Bio import SeqIO

    wanted = set(ordered_ids)
    seq_by_id: dict[str, str] = {}
    for f in sorted((PROC / seq_dir).glob(glob)):
        for rec in SeqIO.parse(f, "fasta"):
            if rec.id in wanted:
                seq_by_id[rec.id] = str(rec.seq)
        if len(seq_by_id) == len(wanted):
            break
    keep = np.array([s in seq_by_id for s in ordered_ids])
    if not keep.all():
        print(f"WARNING: {(~keep).sum()} ids missing in {seq_dir}, dropping", flush=True)
    ordered_ids = [s for s, k in zip(ordered_ids, keep) if k]
    labels = labels[keep]
    seqs = [seq_by_id[s] for s in ordered_ids]
    print(f"found {len(seqs)} sequences in {seq_dir}", flush=True)
    if args.dry_run:
        print("--dry-run OK", flush=True)
        return

    # 3. one pass, all hidden layers, batch_size=1
    import torch
    from darkmatter.device import inference_device
    from darkmatter.embeddings import load_embedder

    preprocess = (lambda s: s)
    if needs_pre:
        from darkmatter.embeddings.genos_m import _preprocess as preprocess

    print(f"loading {args.model}...", flush=True)
    embedder = load_embedder(args.model)
    if args.fp32:
        embedder.model = embedder.model.float()
        print("cast model to fp32 (stable intermediate layers)", flush=True)
    model, tok = embedder.model, embedder.tokenizer
    device = inference_device(embedder.plan, model)

    per_layer: list[list[np.ndarray]] | None = None
    t0 = time.time()
    for j, s in enumerate(seqs):
        toks = tok([preprocess(s)], padding=True, truncation=True,
                   max_length=embedder.max_tokens, return_tensors="pt")
        toks = {k: v.to(device) for k, v in toks.items()}
        with torch.inference_mode():
            hs = model(**toks, output_hidden_states=True).hidden_states
        mask = toks["attention_mask"].unsqueeze(-1)
        if per_layer is None:
            per_layer = [[] for _ in range(len(hs))]
        for li, h in enumerate(hs):
            m = mask.to(h.dtype)
            pooled = (h * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)
            per_layer[li].append(pooled.float().cpu().numpy()[0])
        del hs
        if (j + 1) % 20 == 0:
            print(f"  {j + 1}/{len(seqs)}  ({time.time() - t0:.0f}s)", flush=True)

    layers = np.stack([np.stack(v) for v in per_layer])  # (L+1, n, H)
    print(f"embedded {len(seqs)} seqs across {layers.shape[0]} layers in {time.time() - t0:.0f}s", flush=True)
    out = PROC / f"{args.model}_layer_diagnostic.npz"
    np.savez_compressed(out, layers=layers, labels=labels, ids=np.array(ordered_ids))
    print(f"saved -> {out}", flush=True)

    rows = [(li, separation_gap(layers[li], labels, False), separation_gap(layers[li], labels, True))
            for li in range(layers.shape[0])]
    best = max(rows, key=lambda r: r[2])
    print("\nlayer | gap(raw) | gap(centered)")
    print("------+----------+--------------")
    for li, g_raw, g_cen in rows:
        print(f"{li:5d} | {g_raw:8.4f} | {g_cen:12.4f}{'   <- best' if (li, g_raw, g_cen) == best else ''}")
    print(f"\n{args.model}: last-layer raw {rows[-1][1]:.4f}  ->  best layer {best[0]}, centered {best[2]:.4f}")


if __name__ == "__main__":
    main()
