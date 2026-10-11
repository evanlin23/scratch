#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
R=/home/user/scratch/cs581/gcmvote/results_p35; : > $R/aln.jsonl; : > $R/trees.jsonl; : > $R/magus.jsonl
for n in SIMHIGH_R35 SIMHIGH_R36; do
  [ -f /opt/work/gcmtrees/reps/$n/results.jsonl ] && cat /opt/work/gcmtrees/reps/$n/results.jsonl >> $R/aln.jsonl
  [ -f /opt/work/gcmvote/reps/$n/results.jsonl ] && cat /opt/work/gcmvote/reps/$n/results.jsonl >> $R/aln.jsonl
  [ -f /opt/work/gcmvote/reps/$n/trees.jsonl ] && cat /opt/work/gcmvote/reps/$n/trees.jsonl >> $R/trees.jsonl
  grep "\"dataset\": \"$n\"" /opt/work/gcmtrees/magus.jsonl | grep '"magus"' >> $R/magus.jsonl
done
true
