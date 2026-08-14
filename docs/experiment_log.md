# Experiment Log

One entry per meaningful run. Format:

```
## YYYY-MM-DD — short title
- Commit: <git sha>
- Config: <relevant config/*.yaml values or diff>
- Result: <what happened, key numbers>
- Next: <what this implies for the next step>
```

Entries start once there's a result worth recording (Phase 2 onward).

## 2026-08-13 — full-panel snapshot differencing
- Commit: c4308e4
- Config: R207 panel (502 genomes, phylum-stratified), dark-at-T0 = no Pfam-35 GA hit, characterised-T1-proxy = Pfam-37 hit in a family new since 35.0
- Result: 1,341,100 proteins scanned; 297,798 dark-at-T0 (22.21%); 51,287 characterised-T1-proxy (3.82%); 4,138 positive-proxy (0.309%) -- 41-83x the 50-100 go/no-go floor
- Next: Report to Track 1 for formal go/no-go sign-off; proxy uses Pfam-37 net-new-family, not true InterPro-latest per P1-D3 -- real number should only be larger

## 2026-08-13 — esm2 embedding separation check (genos-m blocked on this gpu)
- Commit: c4308e4
- Config: esm2_t33_650M, mean-pooled; unambiguous single-Pfam-35-hit proteins, 5/family; pilot=1 genome/6 families/24 seqs, full=30 phylum-diverse genomes/20 families/100 seqs
- Result: pilot: within-family cosine 0.972 vs across-family 0.932, gap 0.040. full: within 0.963 vs across 0.928, gap 0.036 -- consistent, modest but real separation
- Next: Genos-m side still pending Angshuman on his M1 Pro; re-run this comparison once his result lands to check separation is comparable across models

## 2026-08-13 — stage 3 embedding sample (esm2, R207 panel)
- Commit: 2ea3410
- Config: esm2, 64000 sequences: {'dark_negative': 30000, 'characterised_at_t0': 29862, 'positive': 4138}
- Result: embedded in 3592s, shape (64000, 1280)
- Next: feed into the Layer 1 novelty scorer once Track 1 fixes EVT vs density-based

## 2026-08-15 — full-panel nucleotide extraction (genos-m input)
- Commit: 9a3e937
- Config: R207 panel (502 genomes), gtdb_proteins_nt_reps streamed directly from GTDB public mirror (61GB archive), no AWS credentials needed
- Result: 502/502 genomes extracted, 1.4GB total, gene counts verified matching the protein version exactly (e.g. 219/219) -- same genes, DNA modality
- Next: Angshuman can now run Genos-m against real, correctly-matched panel data on his M1 Pro
