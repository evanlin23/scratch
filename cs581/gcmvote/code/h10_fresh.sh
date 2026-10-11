#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h10: ROSE 1000M1 R0-R3, one fresh MAGUS draw each (paper flags) -> rep dir -> linsi/es3/es4/es5 baselines.
# Mirrors cs581/gcmgen/code/fresh.sh. Restartable.
W=/opt/work/h10; G=/home/user/scratch/cs581/gcmgen/code
mkdir -p $W/reps $W/fresh
cd /home/user/scratch/cs581/code
for r in ${REPS:-R0 R1 R2 R3}; do
  name=1000M1_$r
  echo "$name /opt/data/Datasets/ROSE/1000M1/$r/rose.aln.true.fasta 25" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads 4 || continue
  python3 /home/user/scratch/cs581/protcons/code/pc.py rep $W/fresh/${name}_d0 $W/reps/$name
  python3 $G/gg.py run $W/reps/$name linsi 'linsi#es3' 'linsi#es4' 'linsi#es5'
done
echo H10_FRESH_DONE
