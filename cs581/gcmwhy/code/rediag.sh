#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Re-run why.py on every rep whose diagnostics predate the h2_split / h6 sections. Restartable.
BANK=/tmp/claude-0/bank; C=/home/user/scratch/cs581/gcmwhy/code; E=$C/../results/edges
while IFS=$'\t' read -r rep typ bb split; do
  [ -e $BANK/$rep ] || continue
  f=$E/$rep.json
  if [ -f $f ] && grep -q h2_split $f; then continue; fi
  BBARG=""; [ "$bb" != "-" ] && BBARG="--bb $BANK/$rep/$bb"
  python3 $C/why.py $BANK/$rep --type $typ $BBARG || echo "FAILED $rep"
done < $C/reps.tsv
echo REDIAG_DONE
