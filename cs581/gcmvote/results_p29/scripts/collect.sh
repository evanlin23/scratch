#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# copy result rows of the p29 replicates into results_p29/
R=$(dirname $(readlink -f $0))/..; : > $R/aln.jsonl; : > $R/trees.jsonl
for n in SIMHIGH_R29 SIMHIGH_R30; do
  cat /opt/work/gcmtrees/reps/$n/results.jsonl /opt/work/gcmvote/reps/$n/results.jsonl >> $R/aln.jsonl 2>/dev/null
  cat /opt/work/gcmtrees/reps/$n/trees.jsonl >> $R/trees.jsonl 2>/dev/null
done
cp /opt/work/gcmtrees/magus.jsonl $R/magus.jsonl
true
