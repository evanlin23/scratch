#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h5: one fresh MAGUS draw (paper flags) per held-out dataset -> rep dir -> gg.py baselines.
# Restartable: bbtool_bench, pc.py rep and gg.py run all skip finished work.
#   bash run_h5.sh NAME TRUE_ALIGNMENT
W=/opt/work/h5; G=/home/user/scratch/cs581/gcmgen/code
name=$1; src=$2
[ -e /opt/work/h5/skip_$name ] && exit 0
mkdir -p $W/reps
cd /home/user/scratch/cs581/code
echo "$name $src 25" > $W/job_$name.txt
python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads 4 || exit 1
python3 /home/user/scratch/cs581/protcons/code/pc.py rep $W/fresh/${name}_d0 $W/reps/$name || exit 1
python3 $G/gg.py run $W/reps/$name linsi 'linsi#es3' 'linsi#es4' 'linsi#es5' || exit 1
echo "BASE_DONE $name"
