#!/bin/bash
# Resilient wrapper for the full-panel 36->37 leakage-tightened boundary scan.
# Unlike the intergenic-control retry wrapper, this job is internally
# checkpointed at batch granularity (label_panel_proteins_36_37.py skips
# every genome already in the output CSV), so a retry must NOT delete
# progress -- it just re-invokes the same command, which picks up wherever
# the last attempt stopped. Loops until the DONE marker appears.
cd /e/dark_matter
DONE_MARKER="data/processed/gtdb_R207/panel_protein_labels_36_37.DONE"
for attempt in $(seq 1 30); do
  echo "=== attempt $attempt/30 at $(date) ===" >> label_36_37_full.log
  uv run python scripts/label_panel_proteins_36_37.py --batch-size 25 --cpus 16 >> label_36_37_full.log 2>&1
  if [ -f "$DONE_MARKER" ]; then
    echo "=== SUCCESS on attempt $attempt ===" >> label_36_37_full.log
    exit 0
  fi
  echo "=== attempt $attempt ended without DONE marker, backing off 30s ===" >> label_36_37_full.log
  sleep 30
done
echo "=== ALL ATTEMPTS FAILED ===" >> label_36_37_full.log
