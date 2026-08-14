#!/bin/bash
cd /e/dark_matter
echo "=== R232 aa-reps (123GB) starting $(date) ==="

attempt=0
max_attempts=200
until uv run python scripts/stream_gtdb_to_s3.py --release R232 --bucket darkmatter-gtdb-067620369122; do
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
