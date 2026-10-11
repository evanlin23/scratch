#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Restartable: one fresh MAGUS draw (paper flags) on RNASim 1000 R0, rep dir, gg baselines + vote variants in 4 lanes.
C=/home/user/scratch/cs581/code; S=/home/user/scratch/cs581; W=/opt/work/gcmvote; name=RNASim
REP=$W/reps/$name; mkdir -p $W/reps
cd $C
if [ ! -d $REP/sets/s0 ]; then
  echo "$name /opt/data/Datasets/RNASim/1000/R0/true_align.txt 25" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/draws --draws 0 --tools '' --e2e '' --union '' --threads 4 || exit 1
  python3 $S/protcons/code/pc.py rep $W/draws/${name}_d0 $REP || exit 1
fi
python3 $S/gcmgen/code/gg.py run $REP linsi 'linsi#es3' 'linsi#es4' 'linsi#es5' > $W/lane_gg.log 2>&1 &
python3 $S/gcmvote/code/run.py $REP magus hard hard-bb > $W/lane_v1.log 2>&1 &
python3 $S/gcmvote/code/run.py $REP es4 soft soft-bb > $W/lane_v2.log 2>&1 &
python3 $S/gcmvote/code/run.py $REP soft2 soft4 > $W/lane_v3.log 2>&1 &
wait
echo H5C_DONE
