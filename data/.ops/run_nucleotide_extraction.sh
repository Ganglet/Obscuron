#!/bin/bash
cd /e/dark_matter
echo "=== R207 nucleotide panel extraction starting $(date) ==="

attempt=0
max_attempts=200
until uv run python scripts/extract_panel_nucleotides.py --release R207; do
    attempt=$((attempt + 1))
    echo "=== python process crashed (attempt $attempt/$max_attempts), $(date) ==="
    if [ "$attempt" -ge "$max_attempts" ]; then
        echo "=== giving up after $max_attempts crashes ==="
        exit 1
    fi
    sleep 60
    echo "=== retrying $(date) ==="
done

echo "=== ALL DONE $(date) ==="
