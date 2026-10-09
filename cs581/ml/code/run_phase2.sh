#!/bin/bash
# Phase 2: wait for phase 1, then IQ-TREE --fast everywhere and the cheap Park-data baselines.
cd "$(dirname "$0")/.."
while pgrep -f "code/run_fast.sh" > /dev/null; do sleep 30; done
python3 code/runtrees.py results/baseline.jsonl --datasets 1000M2 1000M3 1000L1 RNASim --methods iqtree_fast
python3 code/runtrees.py results/baseline.jsonl --datasets Park1000M1HF --reps 0 1 2 3 4 --alns true_align --methods fasttree iqtree_fast
python3 code/runtrees.py results/baseline.jsonl --datasets ParkRNASim1000 --reps 1 2 3 4 5 --alns true_align --methods iqtree_fast
