#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Copy a finished dataset's rows into results_h5: NAME.results.jsonl (gg.py baselines), NAME.vote.results.jsonl
# (run.py vote variants), NAME.<variant>.model.json (vote.py fit, support/exposure cutoffs), NAME.calib.json.
W=/opt/work/h5; D=$(dirname "$0"); name=$1
cp $W/reps/$name/results.jsonl $D/$name.results.jsonl
cp $W/vreps/$name/results.jsonl $D/$name.vote.results.jsonl
cp $W/vreps/$name/calib.json $D/$name.calib.json
for m in $W/vreps/$name/vote/*/model.json; do t=$(basename $(dirname $m)); cp $m $D/$name.$t.model.json; done
grep "\"$name\"" $W/magus.jsonl > $D/$name.magus.jsonl
