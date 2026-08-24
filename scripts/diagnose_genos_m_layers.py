#!/usr/bin/env python
"""Genos-m representation diagnostic: why does last-layer mean-pool collapse?

Re-embeds the SAME 100-seq matched set from compare_embeddings_full.py's
`genos-m_separation_full.json` (reuses its family labels -> NO 13-min Pfam
rescan), capturing ALL hidden layers in one pass, then reports the
separation gap per layer x {raw, mean-centered}. One Genos-m pass; the
per-layer analysis is instant and offline afterwards (raw vectors saved).

Settles whether Genos-m is genuinely weaker than ESM-2 (last-layer gap
0.036) or just badly pooled/anisotropic -- the P1-D13 deferred layer
ablation, triggered because the last-layer separation came back weak
(0.0042, cosines pinned at ~0.99).

    uv run python scripts/diagnose_genos_m_layers.py --dry-run   # validate plumbing, no model
    uv run python scripts/diagnose_genos_m_layers.py             # the ~38-min run (quiet the machine!)

Genos-m is 9.4GB on 16GB unified memory -> quit other apps first or it thrashes swap.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

PROC = Path(__file__).resolve().parents[1] / "data" / "processed" / "gtdb_R207"


def separation_gap(vectors: np.ndarray, labels: np.ndarray, center: bool) -> tuple[float, float, float]:
    """within/across mean cosine and their gap; center=True removes the shared
    mean direction first (the standard anisotropy fix)."""
    v = vectors - vectors.mean(axis=0, keepdims=True) if center else vectors
    v = v / np.clip(np.linalg.norm(v, axis=1, keepdims=True), 1e-12, None)
    sim = v @ v.T
    same = labels[:, None] == labels[None, :]
    off = ~np.eye(len(labels), dtype=bool)
    within = float(sim[same & off].mean())
    across = float(sim[~same & off].mean())
    return within, across, within - across


def load_labeled_subset() -> tuple[list[str], np.ndarray]:
    """Exact 100 ids + family labels from the prior compare run (no rescan)."""
    prior = json.loads((PROC / "genos-m_separation_full.json").read_text())
    families = prior["families"]  # {pfam_acc: [seq_ids]}
    ordered_ids = [sid for ids in families.values() for sid in ids]
    labels = np.array([fam for fam, ids in families.items() for _ in ids])
    print(f"{len(ordered_ids)} sequences, {len(families)} families (reused from prior run)", flush=True)
    return ordered_ids, labels


def load_nt_sequences(wanted_ids: list[str]) -> dict[str, str]:
    from Bio import SeqIO

    wanted = set(wanted_ids)
    seq_by_id: dict[str, str] = {}
    nt_dir = PROC / "panel_nucleotides"
    t0 = time.time()
    for fna in sorted(nt_dir.glob("*.fna")):
        for rec in SeqIO.parse(fna, "fasta"):
            if rec.id in wanted:
                seq_by_id[rec.id] = str(rec.seq)
        if len(seq_by_id) == len(wanted):
            break
    print(f"found {len(seq_by_id)}/{len(wanted)} nt sequences in {time.time() - t0:.0f}s", flush=True)
    return seq_by_id


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true",
                    help="validate label + sequence loading, then exit before loading the model")
    args = ap.parse_args()

    ordered_ids, labels = load_labeled_subset()
    seq_by_id = load_nt_sequences(ordered_ids)

    missing = [s for s in ordered_ids if s not in seq_by_id]
    if missing:
        print(f"WARNING: {len(missing)} ids not found in panel_nucleotides, dropping", flush=True)
        keep = np.array([s in seq_by_id for s in ordered_ids])
        ordered_ids = [s for s, k in zip(ordered_ids, keep) if k]
        labels = labels[keep]
    seqs = [seq_by_id[s] for s in ordered_ids]
    lengths = np.array([len(s) for s in seqs])
    print(f"sequence lengths (nt): min {lengths.min()}, median {int(np.median(lengths))}, "
          f"max {lengths.max()}  (Genos-m truncates at 8192)", flush=True)

    if args.dry_run:
        print("\n--dry-run OK: labels + sequences load cleanly. "
              "Re-run without --dry-run for the full layer sweep.", flush=True)
        return

    # --- one Genos-m pass, all hidden layers, batch_size=1 to cap peak memory ---
    import torch
    from darkmatter.device import inference_device
    from darkmatter.embeddings import load_embedder
    from darkmatter.embeddings.genos_m import _preprocess

    print("loading genos-m...", flush=True)
    embedder = load_embedder("genos-m")
    model, tok = embedder.model, embedder.tokenizer
    device = inference_device(embedder.plan, model)

    per_layer: list[list[np.ndarray]] | None = None
    t0 = time.time()
    for j, s in enumerate(seqs):
        toks = tok([_preprocess(s)], padding=True, truncation=True,
                   max_length=embedder.max_tokens, return_tensors="pt")
        toks = {k: v.to(device) for k, v in toks.items()}
        with torch.inference_mode():
            hidden_states = model(**toks, output_hidden_states=True).hidden_states  # tuple(L+1) (1,T,H)
        mask = toks["attention_mask"].unsqueeze(-1)
        if per_layer is None:
            per_layer = [[] for _ in range(len(hidden_states))]
        for li, h in enumerate(hidden_states):
            m = mask.to(h.dtype)
            pooled = (h * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)  # (1,H)
            per_layer[li].append(pooled.float().cpu().numpy()[0])
        del hidden_states
        if (j + 1) % 10 == 0:
            print(f"  {j + 1}/{len(seqs)}  ({time.time() - t0:.0f}s)", flush=True)

    layers = np.stack([np.stack(v) for v in per_layer])  # (L+1, n_seq, H)
    print(f"embedded {len(seqs)} seqs across {layers.shape[0]} layers in {time.time() - t0:.0f}s", flush=True)

    out = PROC / "genos-m_layer_diagnostic.npz"
    np.savez_compressed(out, layers=layers, labels=labels, ids=np.array(ordered_ids))
    print(f"saved raw per-layer vectors -> {out}", flush=True)

    rows = [(li, *separation_gap(layers[li], labels, center=False)[2:3],
             *separation_gap(layers[li], labels, center=True)[2:3])
            for li in range(layers.shape[0])]
    best = max(rows, key=lambda r: r[2])
    print("\nlayer | gap(raw) | gap(centered)")
    print("------+----------+--------------")
    for li, g_raw, g_cen in rows:
        star = "   <- best" if (li, g_raw, g_cen) == best else ""
        print(f"{li:5d} | {g_raw:8.4f} | {g_cen:12.4f}{star}")
    print(f"\nESM-2 reference (last layer, raw): 0.036")
    print(f"Genos-m last-layer raw: {rows[-1][1]:.4f}  ->  best: layer {best[0]}, "
          f"centered gap {best[2]:.4f}")


if __name__ == "__main__":
    main()
