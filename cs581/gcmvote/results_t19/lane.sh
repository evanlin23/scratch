#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper t19 lane for SIMHIGH_R19: vote.py variants on the gcmtrees MAGUS draw (restartable).
# A separate rep dir under /opt/work/gcmvote/reps links the gcmtrees rep's inputs, because run.py and gg.py both
# append to REP/results.jsonl in different row formats.
src=/opt/work/gcmtrees/reps/SIMHIGH_R19; rep=/opt/work/gcmvote/reps/SIMHIGH_R19
mkdir -p $rep
for f in inputs true.fasta unaligned.fasta magus.json sets; do [ -e $rep/$f ] || ln -s $src/$f $rep/$f; done
G=/home/user/scratch/cs581/gcmvote/code
printf '%s\n' magus es4 hard soft soft2 soft4 hard-bb soft-bb hard+mask | \
  xargs -P ${P:-4} -I{} sh -c "python3 $G/run.py $rep --B 10 '{}' > $rep/log_{}.txt 2>&1"
cat $rep/results.jsonl
