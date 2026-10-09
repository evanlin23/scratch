#!/bin/bash
# Phase 1 (cheap): FastTree baselines everywhere, masking/ensemble pilot alignments, Park RNASim1000.
cd "$(dirname "$0")/.."
python3 code/runtrees.py results/baseline.jsonl --datasets 1000M2 1000M3 1000L1 RNASim --methods fasttree fasttree_pasta
python3 code/runtrees.py results/baseline.jsonl --datasets ParkRNASim1000 --reps 1 2 3 4 5 --alns true_align --methods fasttree
python3 code/runtrees.py results/masking.jsonl --datasets 1000M2 1000M3 1000L1 RNASim \
    --alns gcm_mask50 gcm_mask70 gcm_oracle50 gcm_oracle70 gcm_pasta --methods fasttree
