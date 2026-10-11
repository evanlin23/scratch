#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# One fresh MAGUS draw on 16S.T (paper flags, K=100, vectorized graph builder), then a rep dir as fresh.sh
# builds it. Restartable (bbtool_bench keeps state.json; pc.py rep skips a finished rep).
W=/opt/work/h12; H=/home/user/scratch/cs581/gcmvote/h12
cd /home/user/scratch/cs581/code
python3 $H/memwatch.py $W/magus_mem.json 14400 python3 $H/fastbench.py $W/job.txt $W/magus.jsonl $W/fresh \
  --draws 0 --tools '' --e2e '' --union '' --threads 4 || { echo MAGUS_FAILED; exit 1; }
python3 /home/user/scratch/cs581/protcons/code/pc.py rep $W/fresh/16S.T_R0_d0 $W/reps/16S.T_R0
echo MAGUS_DONE
