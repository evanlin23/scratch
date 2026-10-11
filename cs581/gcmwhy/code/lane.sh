#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# One lane: for every rep in reps.tsv with (line index % NL == L): why.py diagnostics + merge.py variants. Restartable.
L=$1; NL=$2; shift 2; VARS="$*"
BANK=/tmp/claude-0/bank; C=/home/user/scratch/cs581/gcmwhy/code
i=0
while IFS=$'\t' read -r rep typ bb split; do
  if [ $((i % NL)) -eq "$L" ]; then
    BBARG=""; [ "$bb" != "-" ] && BBARG="--bb $BANK/$rep/$bb"
    [ -f $C/../results/edges/$rep.json ] || python3 $C/why.py $BANK/$rep --type $typ $BBARG
    python3 $C/merge.py $BANK/$rep $VARS $BBARG
  fi
  i=$((i+1))
done < $C/reps.tsv
echo LANE_DONE $L
