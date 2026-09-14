#!/usr/bin/env bash
# Resume the ProstT5 dark-query scoring run after any interruption (wifi drop,
# laptop shutdown, OOM kill). Safe to re-run any number of times: completed
# chunks live on disk under data/processed/gtdb_R207/_prostt5_query_chunks/
# and are skipped automatically (see scripts/run_layer2_prostt5.py).
#
# HF_HUB_OFFLINE=1 because both ProstT5 repos are already fully cached
# locally (.cache/huggingface/hub) -- resuming should never need network.
set -euo pipefail
cd "$(dirname "$0")/.."
export HF_HUB_OFFLINE=1
export HF_HUB_DISABLE_XET=1
uv run python scripts/run_layer2_prostt5.py --eval-only --score-queries --batch-size 2 \
    >> prostt5_full_v2.log 2>&1
