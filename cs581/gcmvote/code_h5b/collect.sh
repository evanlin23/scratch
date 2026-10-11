#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Copy one replicate's rows and model stats into results_h5b/: NAME.results.jsonl (vote.py variants, run.py format),
# NAME.gg.jsonl (gg.py baselines), NAME.VARIANT.model.json, NAME.calib.json (calib.py on the magus edges).
name=$1; rep=/opt/work/gcmvote/reps/$name; D=/home/user/scratch/cs581/gcmvote/results_h5b; mkdir -p $D
cp $rep/vote_results.jsonl $D/$name.results.jsonl
cp $rep/results.jsonl $D/$name.gg.jsonl
for d in $rep/vote/*/; do v=$(basename $d); [ -f $d/model.json ] && cp $d/model.json $D/$name.$v.model.json; done
[ -f $rep/calib.json ] && cp $rep/calib.json $D/$name.calib.json
true
