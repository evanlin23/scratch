#!/bin/bash
# Copy per-rep result rows and control rows into the repo.
R=/home/user/scratch/cs581/gcmclust/results
for d in /opt/work/gcmclust/reps/*/; do n=$(basename $d)
  [ -f $d/results.jsonl ] && cp $d/results.jsonl $R/$n.results.jsonl
  [ -f $d/ctrl.json ] && cp $d/ctrl.json $R/$n.ctrl.json
done
[ -f /opt/work/gcmclust/fresh_magus.jsonl ] && cp /opt/work/gcmclust/fresh_magus.jsonl $R/fresh_magus.jsonl
true
