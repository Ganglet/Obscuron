# Phase 1 Track 1 — Retrospective Benchmark & Scope Gate

**Phase:** 1 — Benchmark Construction and Embedding Pipeline (time-boxed, weeks 1–3)
**Track:** 1 — Methodology, Design & Analysis
**Status:** In progress
**Weeks:** 1–3
**Branch:** `phase-1-track-1`

---

## Objective

Turn the blueprint's prose into a defensible, executable specification and pass the go/no-go gate before any modelling begins. Track 1 owns the *what and why*; Track 2 implements the ingestion and snapshot-differencing pipeline against this spec. Nothing is scored yet — the output of this phase is a trusted benchmark definition plus a data-backed decision on whether the primary method is viable.

Five deliverables, in dependency order:

1. Temporal snapshot boundary (Track 1 defines → Track 2 differences).
2. Embedding-extraction protocol (Track 1 specifies → Track 2 runs).
3. Positive-label size estimate + **go/no-go gate** (blocks Phase 2).
4. Systematic literature search (novelty positioning).
5. Documentation & reproducibility standards.

Blueprint references: §2.5 (retrospective validation), §6 Layer 1, §7 (recommended approach, metrics, go/no-go), §8 (work division). Decisions: D2, D3, D4, D5, D6, D7, D10.

---

## 1. Temporal snapshot boundary

**Status:** DECIDED — see P1-D2.

**T₀ = GTDB R207 (Apr 2022) + Pfam 35.0 (Nov 2021); T₁ = GTDB R232 (Apr 2026) + InterPro-latest.**

The choice is driven by **embedding temporal leakage**, not just gap length. Retrospective validity assumes the novelty score reflects only T₀-era knowledge, but the pretrained embedding models postdate T₀, so the embedding space leaks post-T₀ characterisation. The boundary is set so that **ESM-2's UniRef50 training (~2021) aligns with T₀** → the ESM-2 (protein) arm is near leakage-free, while Genos-m (2026) is high-leakage → the ESM-2-vs-Genos-m comparison doubles as a **leakage-sensitivity probe**. Full rationale + rejected alternative (a more recent T₀) in P1-D2.

- [x] Candidate pairs considered; R207→R232 fixed (matches Track 2's fetched snapshots).
- [x] Positive-label count projected (proxy) — §3 / P1-D5.
- [x] T₀/T₁ fixed and recorded (P1-D2).
- [x] Operational definition of "dark at T₀" fixed (§3 / P1-D3).
- [ ] Verify each embedding model's exact training cutoff (feeds the leakage argument).
- [ ] Add one shorter, recent robustness gap (secondary).

**Cross-track hand-off:** Track 2 differences R207-reps against Pfam-35 (dark set) and against InterPro-latest (characterised set); tool = `hmmscan` at GA thresholds.

---

## 2. Embedding-extraction protocol

**Status:** pending

Specify how a sequence becomes a vector — identically for Genos-m and ESM-2 so the comparison (D3) is fair.

**Design note — cutoffs verified; the comparison is a leakage probe (P1-D2, P1-D7).** ESM-2 = UniRef50 **2021_04** (~late 2021) ≈ T₀ → **leakage-clean**, so it carries the retrospective headline. Genos-m = pretrained on **GTDB R220 (24 Apr 2024)** — inside the T₀→T₁ window and the same DB → contaminated; handled by restricting Genos-m positives to **post-R220 (Apr 2024→T₁)** for a clean number, plus the full-window number to quantify contamination. Genos-m leak is self-supervised (weaker) and conservative in direction. Full handling in P1-D7.

To specify:
- [ ] **Model selection & versions** — Genos-m checkpoint (verify against source preprint); ESM-2 size — smoke test uses `esm2_t30_150M`; use a larger variant (`t33_650M`+) for the actual comparison (D9). Record each model's **training-data cutoff**.
- [ ] **Input unit** — nucleotide contig/gene (Genos-m) vs translated protein (ESM-2). Note the asymmetry: Genos-m can see genomic context, ESM-2 sees only the translated ORF. Record the translation/ORF-calling step for the ESM-2 arm.
- [ ] **Layer / pooling strategy** — which hidden layer; mean-pool vs [CLS]/EOS vs per-residue→pooled. Backends currently mean-pool `last_hidden_state`. Keep identical across models where architecture allows.
- [ ] **Batching & precision** — batch size, max sequence length / chunking (Genos-m cap 8192 single-nt tokens ≈ 8kb → long genes need chunking), fp16 on MPS/CUDA.
- [ ] **Normalisation** — L2-normalise embeddings before distance computation? (affects kNN/EVT downstream). Decide and justify.

**Cross-track hand-off:** Track 2 implements extraction to this spec and returns embedding matrices for the labelled subset so Track 1 can compare category separation.

---

## 3. Positive-label size estimate + go/no-go gate

**Status:** PROXY PASSED — gate = **GO** (P1-D3, P1-D5). Definitive count deferred to the Week-2 hmmscan.

**Operational definition (P1-D3):** sequences = GTDB rep-proteins; dark-at-T₀ = no Pfam-35 hit at GA threshold; characterised-by-T₁ = now hits a family new since T₀ (Pfam-37/InterPro-latest); positive = dark ∧ characterised; still-dark = *unlabelled*, not a true negative; GTDB taxonomy is the sequence source, not the function signal.

**Proxy estimate (P1-D5):** Pfam 35.0 = 19,632 families → Pfam 37.0 = 21,979 → **+2,347 net new families**. Each is a dark→characterised cluster → positive population **~10⁴–10⁵**, i.e. ≫ the 50–100 floor (D6) by 2–3 orders. **Availability is not the binding constraint — local embedding compute is (→ subset, §5 / P1-D4).**

- [x] Projected positive count (proxy) — passes with huge margin.
- [x] Compared against floor (D6) — GO.
- [ ] Definitive count: `hmmscan(R207-reps, Pfam-35)` → dark; `hmmscan(dark, InterPro-latest)` → characterised (Week-2, Track 2).
- **Selection-bias note (mandatory, D2):** characterised-since sequences are enriched for the *near-known*, so the positive set skews to moderate novelty, not maximal. Precision@K therefore rewards *real-but-tractable* novelty, not the most distant sequences — state this explicitly in the paper. (Track 1 interpretation deliverable.)

---

## 4. Systematic literature search

**Status:** DONE (first pass) — see P1-D11. Novelty claim narrowed and now defensible; re-run before submission.

Searched web + arXiv + PubMed + bioRxiv for {open-set recognition, novelty/anomaly detection, EVT calibration} × {protein/genomic FM embeddings, metagenomics, microbial dark matter}.

**Closest prior art**
- **HiFi-NN** (iScience 2025) — "annotating the microbial dark matter," ESM-2 650M embeddings. **Annotates** (EC numbers, closed-set); calibration is heuristic kNN-softmax + hard 0.38 cutoff (not conformal/EVT); no open-set rejection; temporal split is a minor benchmark. Its ESM-2 650M choice validates our protein arm.
- **DeepVirus / viral-dark-matter hierarchical DL** (bioRxiv 2025) — open-set + protein FM + genome context + novel-group detection, but **viral lineage** (not microbial function) and hypothesis-testing (not EVT).
- Neighbourhood: protein-embedding anomaly detection (NAR GB 2024); EVM/OpenMax (EVT open-set, general ML); FUGAsseM/CAFA (temporal-holdout function prediction); protein-embedding uncertainty (Nat Methods 2026).

**Refined novelty claim** (replaces the blueprint's over-broad phrasing): *first to cast microbial functional dark-matter prioritisation as formally-calibrated open-set novelty detection (EVT novelty score with a checked error rate, not annotation), validated by a leakage-controlled retrospective benchmark contrasting a leakage-clean protein FM (ESM-2) vs a genomic FM (Genos-m).*

Differentiators: novelty-not-annotation (vs HiFi-NN) · formal EVT calibration (vs HiFi-NN cutoff / DeepVirus hypothesis test) · microbial functional (vs DeepVirus viral) · leakage-controlled retrospective (vs plain temporal holdout). **Caveat to honour: retrospective validation is not itself novel — the leakage control (P1-D7) is.**

- [x] Search run; closest prior art + differentiators recorded (P1-D11).
- [x] Novelty claim narrowed to a defensible statement.
- [ ] Re-run immediately before submission (Phase 4).

---

## 5. Documentation & reproducibility standards

**Status:** in progress.

**Provenance & subsetting standard (P1-D4):** raw reproducible public data (GTDB FASTAs) is **never warehoused** — the pipeline streams needed sequences from the public mirror and persists only *derived* artifacts (embeddings, dark/characterised labels, manifest with source URLs + checksums). The benchmark is a **principled stratified taxonomic subset** (required for local compute — the full set is ~10⁸ proteins; sampling design documented to avoid taxonomic bias). Defer R232 proteins (application-time) and `nt_reps`/`genomes`; lifecycle-expiry on any S3 staging.

- [x] Repo structure, README, LICENSE, `.gitignore`, `ACKNOWLEDGEMENTS.md`, decisions log.
- [x] Provenance + subsetting standard set (P1-D4).
- [ ] Subset spec — how many species/phyla to sample (clears go/no-go with margin, embeddable on M1 Pro). **Next Track 1 task.**
- [ ] Dataset manifest schema — per sequence: id, source snapshot + URL, dark-at-T₀ flag, characterised-by-T₁ flag + label, provenance/checksum.
- [ ] Fixed seeds and version-controlled configs (`config/`) for every run.
- [ ] Branch convention: `phase-N-track-N`; merge to `main` at phase boundaries.

---

## 6. Benchmark subset spec (P1-D6)

The full R207 rep set is ~10⁸ proteins — unembeddable locally and not to be warehoused (P1-D4). The benchmark is a principled subset, built in two stages so genomes are sampled *broadly* (cheap) but only a *budgeted* slice is embedded (the Genos-m bottleneck).

### Stage 1 — genome panel (stratified by phylum)

GTDB R207 species-reps, grouped by GTDB phylum (`taxonomy.csv`). Per phylum:

```
n_phylum = clip( round(f · N_phylum),  floor=2,  cap=30 )      # f ≈ 2%, tuned so Σ ≈ G
G ≈ 500 genomes
```

Scheme = **proportional-with-floor-cap**. Rationale: `floor=2` guarantees every phylum appears — including the rare candidate phyla (CPR/DPANN) where dark matter concentrates, which pure-proportional would bury under Proteobacteria; `cap=30` stops giant phyla dominating. This is the choice that keeps novelty results from being a Proteobacteria artifact — state it in the paper. Fix and record the sampling seed.

### Stage 2 — label the panel (CPU; Track 2)

`hmmsearch` panel proteins against **Pfam-35 profiles** (fast direction) and **InterPro-latest**:
- **dark-at-T₀** = no Pfam-35 hit at GA threshold
- **positive** = dark-at-T₀ ∧ hits an InterPro entry that *postdates T₀* (new-family condition — excludes threshold-artifact false positives where a T₀-era family was merely missed)
- **characterised-at-T₀** = has a Pfam-35 hit → the known-boundary reference set

### Stage 3 — embedding budget (GPU; M1 Pro)

Target **~50–80k** embeddings (Genos-m overnight; ESM-2 faster):
- **all positives** — rare, never subsampled;
- capped **dark-then-still-dark** negatives (~20–40k) — the query universe for Precision@K;
- stratified **characterised-at-T₀** reference (~20–40k, per-phylum cap) — the "known boundary" the novelty score measures distance from.

### Pilot first (de-risk the rate)

Before the full run, label a **50-genome pilot** to *measure* the real dark-fraction and positive-rate on GTDB microbial data (the ~20–40% dark / ~0.5–1% positive-rate figures are UniProt-wide extrapolations). Set G from the measured rate. Expected positives at G=500: ~1,000 (pessimistic 0.3%) to ~2,800 (~0.8%) ≫ the 100 floor — robust to a 3× rate miss.

### Reuse for the held-out-family AUROC

The held-out-family sanity test (blueprint metric 2) reuses this same panel — hundreds of genomes yield thousands of families, so a known family can be withheld from the reference and its members checked for high novelty.

### Cross-track hand-off (Track 2)

Implement against `taxonomy.csv`:
1. Group species-reps by phylum; allocate `n_phylum` (formula above); sample with a fixed seed → panel manifest (genome IDs + accessions + source URLs).
2. Stream each panel genome's proteins from the public GTDB mirror (P1-D4 — no warehousing); `hmmsearch` vs Pfam-35 + InterPro-latest → per-protein labels.
3. Run the 50-genome pilot first; report measured dark-fraction + positive-rate to Track 1 before the full panel.
4. Emit the labelled manifest (id, phylum, source URL, dark/characterised/positive, checksum) — the committed artifact; sequences and embeddings stay gitignored.

Track 1 owns the parameters (f, floor, cap, G, embedding budget, seed) and the allocation scheme; Track 2 owns implementation and execution.

## Decisions logged

- [x] Snapshot pair T₀/T₁ + leakage rationale — **P1-D2**.
- [x] Operational definition of dark / characterised — **P1-D3**.
- [x] Data provenance + subsetting standard — **P1-D4**.
- [x] Go/no-go proxy outcome (GO) — **P1-D5**.
- [x] Benchmark subset spec (stratified panel + budget) — **P1-D6** (§6).
- [x] Model training cutoffs verified + per-arm leakage handling — **P1-D7**.
- [x] Literature search + narrowed novelty claim — **P1-D11** (§4).
- [ ] Embedding layer/pooling/normalisation (§2) — pending full spec.

## Next Track 1 tasks

1. Finalise embedding-extraction spec (§2) — layer/pooling/normalisation, ORF-calling for the ESM-2 arm.
2. Hand the subset spec (§6) to Track 2; review the 50-genome pilot's measured rates before the full run.
3. Re-run the literature search before submission (§4).
