#!/bin/bash
# Runs all replicates (4 in parallel, CAMUS single-threaded); restartable.
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=${OUT:-/opt/runs/camus/res}
jobs() {
  for r in $(seq -w 20 39); do echo "n25 $r g_500"; done
  for r in $(seq -w 0 19); do echo "n15 $r iqtree_500"; done
  for r in $(seq -w 0 19); do echo "n25 $r iqtree_500"; done
  for r in $(seq -w 0 19); do echo "n50 $r g_500"; done
}
jobs | xargs -P 4 -L 1 sh -c 'python3 -I '"$HERE"'/run_rep.py $0 $1 $2 '"$OUT"' >> '"$OUT"'/log_$0_$2_$1.txt 2>&1'
echo ALLDONE
