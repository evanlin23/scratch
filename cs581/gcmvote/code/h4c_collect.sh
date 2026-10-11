#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Copy h4c rows into results_h4c/: NAME.results.jsonl (vote.py rows, rep renamed from NAME_vote to NAME),
# NAME.gg.results.jsonl (gg.py baselines), NAME.VARIANT.model.json, and the fresh MAGUS draw rows.
W=/opt/work/gcmvote; O=/home/user/scratch/cs581/gcmvote/results_h4c; mkdir -p $O
for name in "$@"; do
  sed "s/\"rep\": \"${name}_vote\"/\"rep\": \"$name\"/" $W/reps/${name}_vote/results.jsonl > $O/$name.results.jsonl
  cp $W/reps/$name/results.jsonl $O/$name.gg.results.jsonl
  for d in $W/reps/${name}_vote/vote/*/; do cp $d/model.json $O/$name.$(basename $d).model.json; done
done
cp $W/magus.jsonl $O/fresh_magus.jsonl
