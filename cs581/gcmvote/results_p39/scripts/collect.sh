#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# collect.sh NAME... : results_p39/{aln,trees,magus}.jsonl from the replicates' rows (results_t20 formats)
O=/home/user/scratch/cs581/gcmvote/results_p39; : > $O/aln.jsonl; : > $O/trees.jsonl; : > $O/magus.jsonl
for n in "$@"; do
  cat /opt/work/gcmtrees/reps/$n/results.jsonl /opt/work/gcmvote/reps/$n/results.jsonl >> $O/aln.jsonl
  cat /opt/work/gcmtrees/reps/$n/trees.jsonl >> $O/trees.jsonl 2>/dev/null
  grep "\"dataset\": \"$n\"" /opt/work/gcmtrees/magus.jsonl >> $O/magus.jsonl
done
