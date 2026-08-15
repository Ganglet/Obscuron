# Obscuron — Problems & Decisions Log

Running log of every non-obvious design decision and every problem encountered during the project. Update this file immediately when a new issue arises or a design choice is made. Entries are numbered sequentially; phase is noted in brackets. `D#` = pre-development / blueprint decisions; `P#-D#` = decisions made during phase #.

---

## Pre-Development Decisions (Blueprint, July 2026)

### D1 — Layer 1 (calibrated novelty scoring) as the methodological core [Blueprint §7]
Five candidate architectures were scoped (§6). Layer 1 — calibrated novelty scoring over genomic embeddings — is the core, with Layer 3 (immune-inspired self/non-self) as a supporting theoretical narrative that reuses the Layer 1 embeddings at minimal extra cost. A scoped extension from Layer 2 (sequence–structure fusion) or Layer 4 (multi-signal convergence) is incorporated in the final phase only if time permits. Reason: Layer 1 has the strongest ratio of methodological novelty to what is implementable in one term on the available hardware. Alternative rejected as *primary*: Layer 2 (full structure prediction) — ESMFold is memory-constrained on 8GB VRAM (see D9).

### D2 — Retrospective (temporal) validation as the core evaluation strategy [Blueprint §2.5]
Dark matter has no contemporaneous ground truth, which blocks not just prediction but *evaluation*. Decision: freeze a reference-database snapshot at an earlier time T₀, score what was dark matter at T₀ using only information available then, and grade those scores against the current, more complete database — where some T₀-dark sequences have since been characterised. This manufactures a genuine, falsifiable test set. Known limitation carried forward (see D-risk in §10): the characterised-since set is *not* a random sample of dark matter — it is enriched for sequences that were already close enough to something known to become tractable, so the most genuinely novel cases may never enter the positive set. Stated explicitly, not treated as resolved.

### D3 — Genos-m primary, ESM-2 as formal fallback *and* comparison [Blueprint §7]
Genos-m (genomic foundation model) is the primary embedding source, but it is an unpublished May-2026 preprint — availability, stability, and reported performance cannot be assumed. Decision: adopt ESM-2 (protein language model) as a formal fallback *and* comparison pipeline, not an optional afterthought. Both embeddings are extracted in Phase 1 on a small labelled subset and their separation of known functional categories compared directly. This de-risks the Genos-m dependency and turns it into a research question — does genome-level context add value over an established protein-level model. If Genos-m underperforms or becomes unavailable at any phase boundary, ESM-2 becomes primary with no change to the surrounding methodology. Alternative rejected: Genos-m only — single point of failure on an unreviewed model.

### D4 — Open-set recognition framing; novelty as a calibrated score, not a class [Blueprint §6, Layer 1]
Dark-matter detection is formulated as open-set recognition: rather than forcing each sequence into a known category, model the boundary of characterised sequence space and assign a calibrated novelty score reflecting distance from it. The concrete scorer — extreme-value-theory (EVT) tail calibration vs density-based scoring (kNN / one-class / isolation forest) — is deferred to a Phase 2 decision, but the *framing* (open-set, ranked score, "none of the above" allowed) is fixed here. The project title commits to *calibrated* novelty, which biases the Phase 2 choice toward EVT as the headline with a density-based baseline.

### D5 — Evaluation metrics fixed before implementation [Blueprint §7]
To prevent post-hoc metric fishing, the evaluation metrics and the value(s) of K are fixed prior to Phase 2 and not adjusted based on results: (1) **Precision@K** — fraction of the top-K novelty-ranked candidates that were characterised between the historical and current snapshots; (2) **AUROC on a held-out-family sanity test** — a known family is temporarily withheld and the method's ability to flag it as anomalous is verified; (3) a **calibration check** comparing predicted novelty scores against observed characterisation outcomes. Complementarity: (2) tests sensitivity to novelty, (1) tests that the surfaced novelty is real, (3) tests that the score is honest.

### D6 — Go/no-go dataset gate [Blueprint §7]
The retrospective positive-label set size is estimated in Phase 1, before any modelling. Floor: a positive set of at least **50–100** characterised sequences is treated as the minimum for a stable AUROC estimate. If the projected count falls below this on the initially selected snapshot pair → widen the interval between snapshots, or broaden references (GTDB together with Pfam/UniProt updates), before Phase 2. If neither yields a viable sample → **Layer 5** (statistical discrimination of coding structure from noise, which requires no positive-label set) is promoted to the primary contribution.

### D7 — GTDB as the primary snapshot pair [Blueprint §2.5, §11]
The Genome Taxonomy Database (GTDB), being regularly versioned, supplies the historical/current snapshot pair for retrospective validation. Its release history gives well-defined temporal boundaries and provenance. Pfam / UniProt updates are held in reserve to broaden the characterised set at the go/no-go gate (D6).

### D8 — Target venues fixed in advance [Blueprint §7]
Primary: **ACM-BCB** (ACM SIGBio flagship; Scopus + Web of Science; equivalent 2027 deadline projected ~Feb–Mar 2027, giving ~3–4 months of manuscript runway after Phase 4). Fallback: **IEEE BIBM** (CCF-B; Scopus; 2027 cycle projected ~July 2027). Dates are projected from prior-year cycles and must be confirmed against the official calls for papers.

### D9 — Compute scope: audit-only on M1 Pro + RTX 4060 (8GB) [Blueprint §9]
No model training; forward inference on pretrained models plus lightweight statistical modelling (kNN, EVT, isolation forest) fits both configurations. The one constraint is full 3-D structure prediction (ESMFold) on 8GB VRAM → restrict to a small high-priority subset with sequence chunking, and prefer the ProstT5 approximation. If full-scale structure prediction is required, rent bounded cloud GPU capacity or request university-cluster allocation — raised early in Phase 1, not after the dataset is scoped.

### D10 — Phase 1 time-boxed at three weeks [Blueprint §7]
Phase 1 (reference-snapshot construction + embedding pipeline + go/no-go gate) is fixed at weeks 1–3. If GTDB snapshot-differencing is not complete by the end of week 3, the contingency is to *narrow scope* to a smaller set of taxa or gene families — the remaining schedule is not renegotiated. A partial, narrower benchmark that keeps Phase 2 on schedule is preferred over a complete benchmark that delays it.

---

## Provenance note

An earlier capstone blueprint (MRI operator-shift conformal prediction, `Capstone_Blueprint_v2.pdf`) was **abandoned** before any implementation. Obscuron replaces it; no code or decisions carry over.

---

## Phase 1 Log

_Phase 1 Track 1 working record: `01_track1_phase1_benchmark_scope.md`. Log decisions as `P1-D#` below as they are made._

### P1-D1 — M1 Pro compute environment verified for ESM-2 on MPS; Genos-m deferred to a supervised test / cloud [Phase 1]
Track 2 (Rayyan) verified the environment on the RTX 4060 but could not test the M1 Pro — it is the Track 1 machine, not his. Ran `uv sync` + `scripts/smoke_test.py --skip-genos-m` on the Track 1 MacBook Pro (M1 Pro, **16GB unified memory**): environment builds cleanly (bitsandbytes correctly skipped on Darwin, torch 2.13.0 MPS), and ESM-2 loads + embeds on MPS → `(3, 640)`. Genos-m (fp16, ~9.4GB, no MPS quantization) was **not** run inline: 16GB unified is tight and the 4060 evidence points to OOM. Apple unified memory is a genuinely different regime from the 4060's hard 8GB VRAM wall, so a fit is not foreclosed — left as a deliberate supervised test (~9GB download + freeze risk on the active machine). ESM-2 confirmed as the practical M1 Pro backend, matching blueprint §7/D3; Genos-m at scale stays the cloud/cluster question (D9). Full detail: `docs/reproducibility.md` → "M1 Pro (Apple Silicon) compute environment".

> **Numbering note:** Track 1's `week-01-benchmarking` branch (unmerged as of this writing, pushed 2026-08-10) independently uses `P1-D2` through `P1-D7` for a different set of decisions (snapshot boundary, dark/characterised definition, subset spec, leakage handling). The entries below continue sequentially from what's actually on `main`/this branch and **will collide with that branch's numbering on merge** — renumber one side then, don't treat these as pre-reconciled.

### P1-D2 — Pfam-37 net-new-family proxy stands in for InterPro-latest [Phase 1, Track 2]
Track 1's operational definition (on the unmerged branch, see note above) specifies InterPro-latest as the "characterised by T1" signal. Full InterProScan is a much heavier lift than anything else in the pipeline — ~6.6GB of data across roughly ten bundled member-database search tools (not just HMMER), realistically needing Docker on this Windows machine. Decision: use a Pfam-37 net-new-family proxy instead — a protein counts as characterised-T1-proxy if it hits a Pfam-37 family that's new since Pfam-35, reusing the same mechanism behind Track 1's own go/no-go proxy estimate. Verified Pfam family counts exactly match Track 1's cited numbers (35.0: 19,632 families; 37.0: 21,979). Limitation stated, not buried: this proxy is narrower than true InterPro-latest — it misses cases where a sequence is newly characterised only via a non-Pfam member database, so the real positive count should only be larger. Alternative rejected for now: full InterProScan via Docker — left as a follow-up if the proxy's fidelity becomes a real concern, not attempted first given the setup cost relative to the time available.

### P1-D3 — Genos-m needs nucleotide sequences, not the protein reps already fetched [Phase 1, Track 2]
Problem, not a decision: every sequence fetched and extracted through Week 2 (`gtdb_proteins_aa_reps`) is amino-acid protein — correct for ESM-2 and the Pfam differencing pipeline, but Genos-m does single-nucleotide tokenization and cannot take protein as input. Protein can't be reverse-derived into its source DNA either (the genetic code is degenerate — multiple codons map to the same amino acid), so this needed a fresh fetch, not a transform of what already existed. Flagged externally by Track 1 before it was caught internally. Fix: built a separate extractor (`src/darkmatter/data/panel_nucleotides.py`) pulling `gtdb_proteins_nt_reps` directly from GTDB's public mirror (not through S3 — no AWS credentials needed, consistent with the provenance standard of streaming what's needed rather than warehousing). Confirmed the archive's gene calls match the protein version exactly per genome (e.g. 219/219), so the existing panel's accession list applies unchanged. Ran to completion: 502/502 genomes, 1.4GB.

### P1-D4 — Full-panel runs need fine-grained checkpointing, not just per-file idempotency [Phase 1, Track 2]
Problem: multi-hour full-panel operations (protein extraction, Pfam labeling, nucleotide extraction) were repeatedly killed mid-run by environment/session restarts outside the pipeline's control — three separate times for the labeling step alone, twice more for the nucleotide extraction (once from a confirmed genuine network/DNS outage, not a code bug). Being resumable only at the level of "skip files already fully written" isn't enough when a single run item (one Pfam scan pass, one archive scan) itself takes tens of minutes to hours. Decision: (1) batch multi-hour operations into small chunks (25 genomes) that checkpoint to disk after every chunk, not just at the end of the whole run; (2) for pure network operations, add OS-level crash resilience — a bash retry-loop wrapper, a logon-triggered resume task, and a periodic watchdog that kills and restarts a hung-but-not-crashed process — mirroring the pattern already built in Week 1 for the S3 transfer. Without this, both the full-panel labeling run and the nucleotide extraction would have had to restart from zero after every interruption instead of losing only the in-flight chunk.
