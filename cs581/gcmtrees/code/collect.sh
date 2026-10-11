#!/bin/bash
# Copy small result rows into the repository.
W=/opt/work/gcmtrees; R=/home/user/scratch/cs581/gcmtrees/results
cat $W/reps/*/results.jsonl > $R/aln.jsonl 2>/dev/null
cat $W/reps/*/trees.jsonl > $R/trees.jsonl 2>/dev/null
cat $W/reps/*/iqtrees.jsonl > $R/iqtrees.jsonl 2>/dev/null
[ -f $W/magus.jsonl ] && cp $W/magus.jsonl $R/magus.jsonl
true
