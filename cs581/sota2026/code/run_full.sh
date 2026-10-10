#!/bin/bash
# Full-dataset benchmark (restartable). Run on an otherwise idle machine.
#   bash run_full.sh TOOLS [JOBS] [TIMEOUT]
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=$HERE/../results/full.jsonl
python3 "$HERE/bench.py" "${2:-$HERE/jobs_full.txt}" "$OUT" --work /opt/runs/sota/full --timeout "${3:-1800}" --tools "$1"
