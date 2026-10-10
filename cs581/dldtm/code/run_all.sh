#!/bin/bash
# Run every (condition, replicate) of the pipeline, 4 at a time, 1 thread each. Restartable.
# Usage: bash run_all.sh OUTROOT MAXSUB METHODS
OUT=${1:-/opt/dldtm_runs}; MAX=${2:-50}; M=${3:-IQfull,FT,IQ,BME,NNJ,PF}
cd "$(dirname "$0")"
for r in R0 R1 R2 R3 R4; do for c in 1000M2 1000L1 1000S3 RNASim1000; do echo "$c $r"; done; done |
  xargs -P 4 -L 1 sh -c 'python3 run_rep.py $0 $1 '"$OUT $MAX $M"' > '"$OUT"'/log_$0_$1_'"$MAX"'.txt 2>&1 || echo FAILED $0 $1'
