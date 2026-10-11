#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Copies SIMHIGH_R6 rows (helper t06) into this directory: aln.jsonl = gg.py rows (gcmtrees rep) then vote run.py rows
# (rows with "B"); trees.jsonl = protbench trees.py rows; magus.jsonl = bbtool_bench MAGUS draw row.
O=$(dirname $0); R=/opt/work/gcmtrees/reps/SIMHIGH_R6; V=/opt/work/gcmvote/reps/SIMHIGH_R6
cat $R/results.jsonl $V/results.jsonl > $O/aln.jsonl
[ -f $R/trees.jsonl ] && cp $R/trees.jsonl $O/trees.jsonl
cp /opt/work/gcmtrees/magus.jsonl $O/magus.jsonl
true
