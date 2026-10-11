#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Merge-only variants on the 16S.T rep, one at a time, each under memwatch (4 h limit, peak RSS).
# Usage: variants.sh RUNNER VARIANT...   RUNNER = gg (gcmgen/code/gg.py run) or a python script taking REP VARIANT.
# Restartable: finished variants are skipped by the runner; a killed variant's stale lock is removed first.
W=/opt/work/h12; REP=$W/reps/16S.T_R0; H=/home/user/scratch/cs581/gcmvote/h12
runner=$1; shift
cd /home/user/scratch/cs581/code
mkdir -p $W/mem
for v in "$@"; do
  s=$(echo "$v" | sed 's/[^A-Za-z0-9._-]/_/g')
  if [ -f $REP/results.jsonl ] && python3 -c "import json,sys; sys.exit(0 if any(json.loads(l)['variant']=='$v' for l in open('$REP/results.jsonl')) else 1)"; then continue; fi
  find $REP/variants -name "*.lock" -delete 2>/dev/null
  if [ "$runner" = gg ]; then cmd="python3 /home/user/scratch/cs581/gcmgen/code/gg.py run $REP $v"; else cmd="python3 $runner $REP $v"; fi
  python3 $H/memwatch.py $W/mem/$s.json 14400 $cmd || echo "VARIANT_FAILED $v"
done
echo VARIANTS_DONE
