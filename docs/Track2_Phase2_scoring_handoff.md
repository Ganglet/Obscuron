# Phase 2 Track 2 — Scoring Module Hand-off (interface contract)

**Authored by Track 1** as the concrete contract for the Phase-2 scorer. Track 2
implements against these signatures / shapes / I-O; Track 1 owns the *what and
why* (design doc: `docs/Track1_phase2_scorer_design.md`, decisions P2-D1..D5).
Frozen hyperparameters live in `config/scorer.yaml` — read them, never hard-code.

**Do the ESM-2 arm end-to-end first** (leakage-clean headline, P1-D7). It does not
wait on the Genos-m nt embed-path fix (P1-D12); Genos-m is the second arm, same
code, once its embeddings exist.

---

## Pipeline

```
Stage-3 embeddings (.npy)  +  protein table (labels)
                 │
                 ▼
   NoveltyScorer.fit(reference)          reference = characterised-at-T0
     ├─ leave-one-out kNN distances      queries   = dark-at-T0
     └─ fit GPD tail at u = P90
                 │
                 ▼
   NoveltyScorer.score(queries) → {distance, pvalue, novelty}
                 │
        ┌────────┴─────────┐
        ▼                  ▼
 evaluate_scorer.py   heldout_family_eval.py
  Precision@K + lift    AUROC dist + calibration curve
   (uses positives)      (controlled ground truth)
```

---

## Input the scorer consumes — the "protein table"

Track 2 assembles ONE table (a DataFrame / parquet) keyed by `protein_id`, aligned
to the embedding matrix, by joining the embedding manifest
(`{arm}_panel_embeddings_manifest.csv`) with the Stage-2 label outputs
(`label_panel_proteins.py`). Required columns:

| column          | type        | meaning                                                            |
|-----------------|-------------|--------------------------------------------------------------------|
| `protein_id`    | str (PK)    | GTDB gene id                                                       |
| `emb_row`       | int         | row index into `{arm}_panel_embeddings.npy`                       |
| `role`          | enum        | `reference` (has Pfam-35 hit) or `query` (dark-at-T0) — from category (P1-D3) |
| `is_positive`   | bool        | dark ∧ characterised-by-T1 (P1-D12); meaningful for `query` rows only |
| `pfam_family`   | str \| null | family id (for `reference`: the held-out unit, P2-D4; null for dark) |
| `phylum`        | str         | GTDB phylum (stratification + matched controls)                    |
| `length_aa`     | int         | for the ≤1024 aa fair-comparison filter + length-matched controls  |

Embeddings are already L2-normalised (P1-D13); if not, normalise on load. The
scorer only ever sees the matrix + this table.

---

## Module API — `src/darkmatter/scoring/`

### `knn.py`
```python
def knn_distance(
    query: np.ndarray,          # (Nq, D) L2-normalised
    reference: np.ndarray,      # (Nr, D) L2-normalised
    k: int = 5,
    reduce: str = "mean",       # mean of the k nearest (P2-D2)
    exclude_self: bool = False, # True => leave-one-out (reference-vs-reference calibration)
) -> np.ndarray:                # (Nq,) score; distance = 1 - cosine = 1 - dot (vectors are unit-norm)
```
`exclude_self=True` drops the zero-distance self-match when `query` rows are also
`reference` rows (the calibration pass). Use an ANN backend (faiss/sklearn) if the
brute-force `Nq×Nr` is too big; result must be identical to brute force.

### `gpd.py`
```python
@dataclass
class GPDTail:
    u: float        # threshold distance (P90 of reference LOO distances)
    zeta_u: float   # empirical P(distance > u) = fraction of reference exceeding u
    xi: float       # GPD shape
    beta: float     # GPD scale

def fit_gpd_tail(ref_distances: np.ndarray, threshold_quantile: float = 0.90) -> GPDTail
    # u = quantile(ref_distances, q); fit scipy.stats.genpareto to (exceedances - u) via MLE

def tail_pvalue(model: GPDTail, distances: np.ndarray) -> np.ndarray:
    # d > u : p = zeta_u * (1 - genpareto.cdf(d - u; c=xi, scale=beta))   (calibrated false-flag rate)
    # d <= u: p from the empirical CDF of ref_distances (below-tail region, not the headline)

def novelty(model: GPDTail, distances: np.ndarray) -> np.ndarray:   # s = 1 - p

def diagnostics(ref_distances: np.ndarray, model: GPDTail) -> dict:
    # mean_residual_life points, QQ vs fitted GPD, Anderson-Darling stat (P2-D5)
```

### `scorer.py`
```python
class NoveltyScorer:
    def __init__(self, config: dict): ...       # parsed config/scorer.yaml

    def fit(self, reference: np.ndarray) -> None:
        # 1. LOO kNN distances over reference (knn_distance(..., exclude_self=True))
        # 2. self._tail = fit_gpd_tail(loo_distances, threshold_quantile)

    def score(self, query: np.ndarray) -> pd.DataFrame:
        # columns: distance, pvalue, novelty   (headline: GPD-calibrated)

    def score_baseline(self, query: np.ndarray) -> np.ndarray:
        # raw kNN distance only, no GPD (the HiFi-NN-style baseline, P2-D1)
```

---

## Scripts

### `scripts/evaluate_scorer.py  --arm esm2 --config config/scorer.yaml`
1. Load matrix + protein table; split `reference` / `query` by `role`.
2. `NoveltyScorer(config).fit(reference)`; `.score(query)` and `.score_baseline(query)`.
3. **Precision@K** at each `K` on the novelty-ranked queries; **lift = P@K / (positives / len(queries))**; report the exact set positive-rate and the full P@K curve. Same for the baseline ranking (the number the headline must beat).
4. Write `results/scorer_{arm}_precisionK.csv`, the ranked query list, GPD params + diagnostics; call `log_experiment(...)`.

### `scripts/heldout_family_eval.py  --arm esm2 --config config/scorer.yaml`
1. Candidate families = `reference` Pfam families with ≥ `min_members`; if > `max_families`, stratified sample by `stratify_by` (seed from config).
2. For each family `F`: remove `F` members from the reference; refit the scorer on reference-minus-`F`; score `F` members (label **novel**) and a matched sample of retained knowns (label **known**, matched on `controls_match`).
3. Aggregate: **AUROC** per family → distribution; pool predicted p-values → **reliability diagram** (calibration, `bins` from config).
4. Write `results/heldout_{arm}_auroc.csv`, calibration bins, and the plots; `log_experiment(...)`.

---

## Acceptance sanity (ESM-2 arm — expected before trusting the numbers)

- Reference LOO distances have a right tail the GPD fits: finite `xi`, `beta`; AD / QQ not catastrophic.
- Held-out-family **AUROC median clearly > 0.5** — the scorer detects genuine novelty (the separation gap was 0.036 in Phase 1, so expect signal but not a landslide).
- Calibration curve tracks the diagonal — the p-value's "checked error rate" is honest (the headline claim).
- Precision@K **lift > 1** and roughly monotone-decreasing in `K`; interpret against the enriched set base rate (P2-D3), not the natural ~1.4%.

Any of these failing is a *finding*, not a bug to hide — log it and raise with Track 1 before touching the frozen config.

---

## Track 2 to confirm / fill

- Exact Stage-3 embedding artifact paths (the config has best-guess `[TODO]`s).
- The exact recorded Phase-1 panel seed → `config/scorer.yaml: seed`.
- Whether `phylum` / `length_aa` / `pfam_family` are already in the manifest or need joining from the label outputs.
- Genos-m arm is **blocked** until the nt embed-path fix (P1-D12) produces `genos-m_panel_embeddings.npy`.
