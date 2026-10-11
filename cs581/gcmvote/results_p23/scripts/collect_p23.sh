#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# copy rows for SIMHIGH_R23/R24 into results_p23
P=/home/user/scratch/cs581/gcmvote/results_p23; G=/opt/work/gcmtrees; V=/opt/work/gcmvote
: > $P/aln.jsonl; : > $P/trees.jsonl
for n in SIMHIGH_R23 SIMHIGH_R24; do
  cat $G/reps/$n/results.jsonl $V/reps/$n/results.jsonl >> $P/aln.jsonl 2>/dev/null
  cat $G/reps/$n/trees.jsonl >> $P/trees.jsonl 2>/dev/null
done
cp $G/magus.jsonl $P/magus.jsonl
true
