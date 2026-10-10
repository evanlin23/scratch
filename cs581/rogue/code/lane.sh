#!/bin/bash
# usage: lane.sh ALIGNERS INSTANCE...   (restartable; rows appended to ../results/rows_<ALIGNERS>.jsonl)
cd "$(dirname "$0")"
A=$1; shift
for i in "$@"; do
  python3 run.py /tmp/claude-0/-home-user-scratch/690a2cdc-90a2-5a40-b34f-c88f48d84cea/scratchpad/inst/$i ../results/rows_$A.jsonl /tmp/claude-0/-home-user-scratch/690a2cdc-90a2-5a40-b34f-c88f48d84cea/scratchpad/work --aligners $A --tree-workers 3 --use-detectors "${USE_DET-ts,pd,hmm}"
done
