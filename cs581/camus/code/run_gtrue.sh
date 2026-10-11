#!/bin/bash
# Control: true gene trees, n25 reps 20-39 (after runall.sh finishes).
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=/opt/runs/camus/res_gtrue; mkdir -p $OUT
while pgrep -f runall.sh >/dev/null; do sleep 30; done
for r in $(seq -w 20 39); do echo "n25 $r g_true"; done | xargs -P 4 -L 1 sh -c 'python3 -I '"$HERE"'/run_rep.py $0 $1 $2 '"$OUT"' --variants default,true_major,tqmc,swap,z3only,z5only,t0.2 >> '"$OUT"'/log_$1.txt 2>&1'
python3 -I $HERE/objective_check.py $OUT n25 g_true $OUT/objective.tsv
echo GTRUEDONE
