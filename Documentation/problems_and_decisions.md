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

<!-- P1-D1 — ... -->
