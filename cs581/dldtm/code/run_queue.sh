#!/bin/bash
# Throttled queue: start each (cond, rep) job when fewer than 4 run_rep.py processes are running.
# Usage: bash run_queue.sh OUTROOT MAXSUB METHODS "REPS" "CONDS"
OUT=$1; MAX=$2; M=$3; REPS=${4:-"R1 R2 R3 R4"}; CONDS=${5:-"1000M2 1000L1 1000S3 RNASim1000"}
cd "$(dirname "$0")"
for r in $REPS; do for c in $CONDS; do
  while [ "$(pgrep -c -f '^python3 run_rep.py')" -ge 4 ]; do sleep 15; done
  (python3 run_rep.py $c $r $OUT $MAX $M > $OUT/log_${c}_${r}_${MAX}.txt 2>&1 || echo FAILED $c $r) &
  sleep 2
done; done
wait
