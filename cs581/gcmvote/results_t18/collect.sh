#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Collect SIMHIGH_R18 rows (helper t18): gg.py rows (gcmtrees rep) + vote.py rows (run.py), trees, MAGUS draw.
D=/home/user/scratch/cs581/gcmvote/results_t18
cat /opt/work/gcmtrees/reps/SIMHIGH_R18/results.jsonl /opt/work/gcmvote/reps/SIMHIGH_R18/results.jsonl > $D/aln.jsonl
cat /opt/work/gcmtrees/reps/SIMHIGH_R18/trees.jsonl /opt/work/gcmvote/reps/SIMHIGH_R18/trees.jsonl > $D/trees.jsonl 2>/dev/null
cp /opt/work/gcmtrees/magus.jsonl $D/magus.jsonl
for v in /opt/work/gcmvote/reps/SIMHIGH_R18/vote/*/; do cp $v/model.json $D/raw/SIMHIGH_R18.$(basename $v).model.json; done
