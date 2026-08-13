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
