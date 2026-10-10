#!/bin/bash
# Copy small per-dataset result files from the work dir into results/.
R=/home/user/scratch/cs581/protcons/results
W=/opt/work/protcons
cp $W/magus.jsonl $R/magus.jsonl 2>/dev/null
for d in $W/reps/*/; do
  n=$(basename $d)
  for k in results.jsonl gate.json magus.json pairdiff.jsonl trees.jsonl; do [ -f $d/$k ] && cp $d/$k $R/$n.$k; done
done
