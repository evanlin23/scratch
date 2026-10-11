#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# collect.sh NAME : append NAME's rows to results_p41/{magus,aln,trees}.jsonl (gg rows, then run.py vote rows)
set -eu
N=$1; O=/home/user/scratch/cs581/gcmvote/results_p41
grep "\"dataset\": \"$N\"" /opt/work/gcmtrees/magus.jsonl >> $O/magus.jsonl
cat /opt/work/gcmtrees/reps/$N/results.jsonl /opt/work/gcmvote/reps/$N/results.jsonl >> $O/aln.jsonl
cat /opt/work/gcmtrees/reps/$N/trees.jsonl >> $O/trees.jsonl
