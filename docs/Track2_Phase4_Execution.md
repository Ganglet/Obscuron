# Full-Scale Genos-m Eval, Optional Extensions & the ProstT5 Infrastructure Saga

**Phase:** 4 — Extension, Analysis & Manuscript
**Owner:** Track 2 (Rayyan)
**Branch:** `phase-4-track-2-extension-manuscript-prep` (P4-D7/D8), `phase-4-track-2-optional-extras` (P4-D9/D10)
**Status:** Complete. All four items closed: the O(N^2) eval-code fix that unblocked full-scale Genos-m, the full 903-family Genos-m cloud run and its use as a 4th convergence axis, and all three of the optional non-blocking items (true-intergenic Layer-5 controls, the standalone benchmark release, and ProstT5 at full scale as a 5th convergence axis). See P4-D7 through P4-D10 in `docs/problems_and_decisions.md` for the full decision-log entries this doc summarizes and expands on. Track 1's Phase-4 design and build work (Layers 1-6, P4-D1..D6) lives in `docs/Track1_phase4_extension_design.md`.

## Objective

Track 1 built all five blueprint layers (P4-D1..D6) and left two items explicitly cloud-gated (full-family Genos-m eval, Genos-m as a Layer-4 convergence axis) plus three further items explicitly optional and non-blocking (ProstT5 at full scale + as a 5th axis, true-intergenic Layer-5 controls, a standalone benchmark release). Track 2's Phase-4 work closes all of these — the cloud-gated pair first (P4-D7/D8), then the three optional items (P4-D9/D10) once the manuscript's Results sections made clear which numbers were still missing.

## What Was Done

### Unblocking full-scale Genos-m (P4-D7)

- **Found the actual blocker**: P2-D11 flagged the full 903-family Genos-m run as blocked because `heldout_auroc`'s per-family loop recomputed a full reference-vs-reference kNN from scratch for every held-out family — an O(F·N²) cost that was fine at the 300-family cap used for the matched comparison, but not at 903+ families.
- **Rewrote the algorithm, not just the constant**: precompute each reference point's top-`m` nearest neighbours once (chunked over query rows so peak memory is `chunk × N`, never the full `N × N` matrix — the same memory discipline the original per-family `knn()` already used, just paid once instead of `F` times), then for each held-out family filter out any precomputed neighbour that falls inside it and take the top-`K` of what's left. `top_m=500` comfortably exceeds every family size in this benchmark.
- **Validated before trusting it**: the raw (uncentered) path reproduces the original per-family implementation bit-for-bit on synthetic reference sets (exact match to 1e-9). The centered path uses a deliberate approximation (global reference mean instead of a per-fold family-excluded mean) whose real-scale error P2-D9 had already measured directly (~0.001-0.003 AUROC) — small synthetic tests showed more divergence (~0.09) than that, but this is expected (tiny N exaggerates the effect of excluding one family from the mean) and not a red flag, since the number that matters is the one already measured at real scale.
- **A design flaw caught by a real memory-pressure incident**: the first version of this optimisation materialised the full `N × N` similarity matrix once — faster, but reintroducing exactly the memory blow-up the original code's own docstring warned against. Local validation against the real 15,360-point reference hit genuine memory pressure (available RAM fell to ~870MB) before this was caught, which motivated reworking to the chunked top-`m` design before ever running it for real, rather than after spending cloud money on a flawed version.

### Full 903-family Genos-m cloud run + 4th convergence axis (P4-D8)

- **Ran on `g5.xlarge`** (A10G, Obscuron AWS account) using the now O(N²) `reembed_eval.py` and the existing (built-but-never-run, since it needed a GPU Genos-m could load) `embed_genos_m_dark_queries.py`.
- **Result: the family-capping inflation effect generalises.** Full 903-family Genos-m AUROC (layer 9) is 0.665 mean / 0.703 median — markedly lower than the 300-largest-family-capped estimate of 0.739/0.795 from P2-D11. This is the same effect P2-D9 already found for ESM-2 (full-903 0.962 vs a capped-family eval reading higher): capping to the largest families is systematically easier (more reference members per fold), so it inflates the apparent score. Read as a general lesson: any family-capped proxy in this project is an upper bound, not the reportable number.
- **Extended `run_convergence.py` to accept a repeatable `--extra-axis` flag** so composition and Genos-m (and later ProstT5, P4-D10) could combine with the two frozen axes in one run.
- **The headline new finding, reported with an explicit caveat rather than folded into the claim**: Genos-m's top-decile novelty has positive-rate lift 1.30× — enriched for eventual characterisation, the only axis in the whole project (Layer 1, Layer 3, or either of the other two Layer-4 axes) that points this direction. Traced this to Genos-m's pretraining overlap with this benchmark's T0→T1 window (P1-D7/P2-D7's leakage lens) rather than reporting it as a clean fourth line of evidence: a gene Genos-m partially memorised during pretraining, later characterised, would register as elevated "novelty" that is really familiarity-with-what-becomes-known.
- **4-way convergent set**: 35 genes, lift 0.94× (near-neutral) — Genos-m's enrichment partially cancels the other three axes' depletion, consistent with the table above rather than contradicting it.
- **Cost**: ~75 minutes of `g5.xlarge` wall time (11 min reference embed + full-family eval + 39 min dark-query embed + convergence merge), instance terminated immediately after — roughly $1.25 at on-demand pricing.

### True-intergenic Layer-5 controls + standalone benchmark release (P4-D9)

- **`scripts/extract_intergenic_controls.py`**: the shuffled/Markov-1 nulls Track 1 built (P4-D5) are model-free substitutes for a real non-coding control. Built the real thing by streaming GTDB's full genome-assembly archive (65GB, R207) and extracting sequence between consecutive gene calls, using coordinates already embedded in the Prodigal-style FASTA headers this project already had locally.
- **Scope decision, documented**: the archive is organised by accession-number path, not front-loaded by panel membership, so covering all 502 panel genomes would mean streaming close to the full 65GB. Streamed a bounded prefix instead (80-genome / 12GB / 30-minute caps, whichever binds first) and used whichever panel genomes were encountered — same documented-tradeoff pattern as P1-D6/P2-D10/P3-D7.
- **Recovered cleanly from a real network failure**: the streaming connection dropped after 2.26GB (`IncompleteRead`, a live interruption, not a bug). Wrapped the whole script in a bash retry loop (up to 5 attempts, P1-D10's resilience pattern) rather than requiring a manual restart; succeeded on the immediate retry.
- **Result** (44 genomes, byte budget was the binding constraint at 12.0GB/1228s, 81,286 intergenic segments ≥60bp): real dark-gene CPBB median 0.0721 (matches P4-D5's number exactly, a consistency check) vs true-intergenic median 0.0308 — markedly higher than the artificial nulls' 0.011, because real intergenic DNA carries its own biological structure. Coding-vs-true-intergenic AUROC = 0.78, lower than the 0.94 against artificial nulls — reported as the more honest number, not a weaker result.
- **`scripts/package_benchmark_release.py` → `benchmark_release/`**: a self-contained package (README, checksummed manifest, genome panel, protein labels) usable without the rest of the repository — ships only derived labels, never raw sequences, per the project's data-provenance standard (P1-D4).

### ProstT5 full-scale + 5th convergence axis (P4-D10)

By far the longest of the four items, entirely from infrastructure fragility on an 8GB local GPU rather than anything about the underlying result.

- **Missing module recovered**: `src/darkmatter/embeddings/prostt5.py`, referenced throughout `scripts/run_layer2_prostt5.py`, had never actually been committed. Rewritten from scratch: `T5EncoderModel` built from the tokenizer/config repo, ProstT5's separately-hosted plain-fp16 `pytorch_model.bin` loaded by hand into it (`transformers` only auto-loads safetensors), rare residues (U/Z/O/B) mapped to X and space-joined per the ProtT5-family tokenizer convention.
- **OOM-on-load fixed properly, not papered over**: the naive build-fp32-then-cast-then-`.to(device)` pattern fragmented the allocator and OOM'd on transfer. A quick `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` attempt made it worse (segfault). Fixed by building the model on `torch.device("meta")` then `to_empty(device=...).to(target_dtype)` — constructs directly on GPU in the target dtype, no CPU fp32 intermediate, no fragmenting per-parameter transfer.
- **5-hour silent run root-caused and fixed**: the dark-query scoring path embedded all 34,138 queries in one unsorted call with no progress output. Random-length batches each pay the padding cost of their longest member, scattered throughout rather than concentrated at the end. Fixed the same way the reference-embedding path already was: sort by length first, chunk, print progress per chunk.
- **Cloud attempt tried and abandoned cleanly**: with the local run's ETA unclear, tried a parallel AWS `g4dn.xlarge`/T4 attempt — hit the same missing-module problem (the local fix hadn't been pushed before the clone). By the time that was fixed, local was already progressing well past where cloud would have started from cold; cloud was terminated unused (~1.6h idle, ~$0.85).
- **Two more OOM crashes, at the same point, cost real lost progress — until real checkpointing was added.** Both a batch-size-8 and a batch-size-2 attempt crashed with CUDA OOM at almost exactly 96.7% through (33,024-33,536 of 34,138 queries, the longest length-sorted sequences finally exceeding available VRAM). Because per-chunk embeddings only lived in memory and were concatenated once at the very end, each crash lost the entire run. Fixed with real resilience rather than a third guess at the batch size: each 256-sequence chunk now saves to disk as it completes and is skipped on restart; an OOM on a chunk triggers `torch.cuda.empty_cache()` and an automatic retry at half the batch size (up to 4 attempts); `scripts/resume_prostt5.sh` sets `HF_HUB_OFFLINE=1` for restarts, since both ProstT5 HF repos are already cached locally and a resume after a wifi drop or shutdown should never need network. The next attempt passed the 96.7% mark cleanly and finished all 34,138 queries in 2,999s.
- **Result, reported honestly rather than promoted**: ProstT5's structural novelty, added as a 5th Layer-4 axis, correlates with the compositional axis at Spearman ρ=+0.64 — by far the largest off-diagonal value in the convergence table (next-highest is 0.32) — and anti-correlates with Genos-m at ρ=-0.43. Its top-decile novelty depletes positives (lift 0.28×), the same direction as composition (0.26×), consistent with the two axes carrying overlapping rather than independent signal. The five-way convergent set shrinks to 20 genes (lift 1.24×) from the four-way set's 35 genes (lift 0.94×) — consistent with a stricter intersection, not with ProstT5 adding independent corroboration. **Decision: report the 5-axis number for completeness but do not promote it over the 3-axis headline result**, the same standing rule P4-D8 applied to Genos-m's 4th-axis inversion.

## Why (Key Decisions)

**Why fix the eval algorithm instead of just running the O(F·N²) version overnight on a bigger cloud box?**
The per-family cost scales with the number of families, so at 903 families (vs the 300-family cap it was written for) the wall-clock cost would have been roughly 3× worse on top of an already-heavy per-family recomputation — throwing more compute at an asymptotically wrong algorithm is the wrong fix when the algorithm itself has an O(N²) redundancy that a one-time precomputation removes for free.

**Why validate the rewrite on synthetic data before trusting it on the real reference?**
A bit-for-bit match on synthetic cluster configurations is a much stronger correctness check than "the final number looks plausible" — it catches an algorithmic bug regardless of whether that bug happens to move the real-scale number by a little or a lot.

**Why report Genos-m's and ProstT5's convergence-axis anomalies instead of quietly dropping them?**
Both this project's own anti-fishing discipline (P2-D5/P3-D5) and its consistent practice (P3-D7's Layer-3 inversion, P1-D7's leakage caveats) treat "report every axis you tried, including the ones that behaved unexpectedly" as a trust-building discipline, not optional transparency. A convergence result that quietly kept only the axes that agreed would be a weaker paper even with a bigger number.

**Why add per-chunk checkpointing only after the second OOM crash, not the first?**
The first crash (batch_size=8) looked like a batch-size problem and was treated as one — halving the batch size is the standard, cheap first fix for a CUDA OOM. The second crash, at almost exactly the same point even at batch_size=2, was the signal that the actual defect was structural (no checkpointing meant every crash cost 100% of the run, not just the failing chunk) rather than a batch-size tuning problem, and that's what got fixed properly.

## Outputs

| Output | Description |
|---|---|
| `scripts/reembed_eval.py` | O(N²) held-out-family AUROC rewrite (P4-D7), committed. |
| `results/genos-m_novelty_scores.csv` | Full 34,138-query Genos-m dark-query novelty scores (P4-D8). |
| `src/darkmatter/embeddings/prostt5.py` | ProstT5 embedder module, written from scratch this phase (P4-D10). |
| `scripts/run_layer2_prostt5.py` | ProstT5 driver, patched with length-sorted chunking, per-chunk checkpointing, and OOM backoff. |
| `scripts/resume_prostt5.sh` | Offline-safe resume wrapper (`HF_HUB_OFFLINE=1`) for interrupted ProstT5 runs. |
| `data/processed/gtdb_R207/prostt5_query_novelty.csv` | Full 34,138-query ProstT5 structural-novelty scores (P4-D10). |
| `data/processed/gtdb_R207/convergence_results.json` | Layer-4 convergence result, now 5 axes. |
| `scripts/extract_intergenic_controls.py` | True-intergenic control extractor (P4-D9). |
| `data/processed/gtdb_R207/intergenic_controls.csv`, `statistical_results_intergenic.json` | True-intergenic control results. |
| `scripts/package_benchmark_release.py`, `benchmark_release/` | Standalone benchmark release package (P4-D9). |

**Next:** all four items are folded into `paper_draft/manuscript_draft.tex` (Layer 4 section, Limitations, Future Work, Discussion, Conclusion — kept out of git per the project's manuscript-drafts policy). Only remaining project work is manuscript polish (user-led). Track 1 declined to re-freeze Layer 3's detector-count default (P3-D8), so the frozen 5,000-detector default stands.
