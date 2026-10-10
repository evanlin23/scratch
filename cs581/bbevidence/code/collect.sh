#!/bin/bash
# Copy the small per-replicate result files from the work dirs into cs581/bbevidence/results/.
R=$(cd "$(dirname "$0")/../results" && pwd)
for d in /home/user/work/reps/*/; do
  n=$(basename $d)
  for k in results diag; do [ -f $d/$k.jsonl ] && cp $d/$k.jsonl $R/$n.$k.jsonl; done
  for t in $d/aligned/*/s*.json; do [ -f "$t" ] && mkdir -p $R/timing && cp $t $R/timing/${n}__$(basename $(dirname $t))__$(basename $t); done
done
