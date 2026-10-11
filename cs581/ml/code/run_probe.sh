#!/bin/bash
cd "$(dirname "$0")/.."
python3 code/runtrees.py results/timing_probe.jsonl --datasets 1000M3 --reps 0 --alns gcm --methods raxmlng_ft raxmlng iqtree_fast --workers 3
