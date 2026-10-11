#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
O=/home/user/scratch/cs581/gcmvote/results_p21; mkdir -p $O; : > $O/aln.jsonl; : > $O/trees.jsonl
for n in SIMHIGH_R21 SIMHIGH_R22; do
  cat /opt/work/gcmtrees/reps/$n/results.jsonl /opt/work/gcmvote/reps/$n/results.jsonl >> $O/aln.jsonl 2>/dev/null
  cat /opt/work/gcmtrees/reps/$n/trees.jsonl >> $O/trees.jsonl 2>/dev/null
done
cp /opt/work/gcmtrees/magus.jsonl $O/magus.jsonl
true
