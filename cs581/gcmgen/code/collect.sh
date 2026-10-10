#!/bin/bash
# Copy per-replicate result files from the work dirs into cs581/gcmgen/results/.
R=/home/user/scratch/cs581/gcmgen/results
for d in /opt/work/gcmgen/reps/*/; do
  n=$(basename $d)
  for k in results agree; do [ -f $d/$k.jsonl ] && cp $d/$k.jsonl $R/$n.$k.jsonl; done
  [ -f $d/magus.json ] && cp $d/magus.json $R/$n.magus.json
  for t in $d/aligned/*/s0.times.json; do [ -f "$t" ] && mkdir -p $R/timing && cp $t $R/timing/${n}__$(basename $(dirname $t)).times.json; done
done
true
