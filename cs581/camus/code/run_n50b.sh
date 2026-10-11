#!/bin/bash
# n50 reps 12-19 with the preregistered variants only (CPU budget), after current n50 jobs end.
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=/opt/runs/camus/res
while pgrep -f "run_rep.py n50" >/dev/null; do sleep 30; done
for r in $(seq 12 19); do echo "n50 $r g_500"; done | xargs -P 4 -L 1 sh -c 'python3 -I '"$HERE"'/run_rep.py $0 $1 $2 '"$OUT"' --variants default,z3only,z5only,tqmc,true_major,wastral,swap >> '"$OUT"'/log_$0_$2_$1.txt 2>&1'
echo N50BDONE
