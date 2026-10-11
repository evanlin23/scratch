#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h5d, RNASim_R1: fresh MAGUS draw (as gcmgen/code/fresh.sh), gg.py baselines, vote.py variants. Restartable.
set -u
S=/home/user/scratch/cs581; C=$S/code; W=/opt/work/gcmgen; V=/opt/work/gcmvote/reps
name=RNASim_R1; rep=$W/reps/$name; vrep=$V/$name
mkdir -p $W $V
cd $C
if [ ! -d $rep/sets/s0 ]; then
  echo "$name /opt/data/Datasets/RNASim/1000/R1/true_align.txt 25" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads 4 || exit 1
  python3 $S/protcons/code/pc.py rep $W/fresh/${name}_d0 $rep || exit 1
fi
# vote replicate: same inputs via symlinks (run.py and gg.py rows go to separate results.jsonl)
if [ ! -d $vrep ]; then
  mkdir -p $vrep.tmp && for f in $rep/*; do ln -s $f $vrep.tmp/; done && mv $vrep.tmp $vrep
fi
# 4 lanes in parallel (gg.py locks per variant; run.py lanes take disjoint variants)
G=$S/gcmgen/code/gg.py; R=$S/gcmvote/code/run.py
( python3 $G run $rep linsi; python3 $G run $rep 'linsi#es3' ) > $W/lane1_$name.log 2>&1 &
( python3 $G run $rep 'linsi#es4'; python3 $G run $rep 'linsi#es5' ) > $W/lane2_$name.log 2>&1 &
( python3 $R $vrep magus es4 hard soft ) > $V/lane3_$name.log 2>&1 &
( python3 $R $vrep soft2 soft4 hard-bb soft-bb ) > $V/lane4_$name.log 2>&1 &
wait
echo PIPELINE_DONE
