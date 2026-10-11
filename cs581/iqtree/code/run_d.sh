#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Per replicate: vote.py merges (magus, es4, hard-bb, scored by gcmvote/code/run.py), then IQ-TREE on
# true/magus/es4/hard-bb, 2 runs at a time. Restartable. Usage: bash run_d.sh SIMHIGH_R7 [...]
set -u
H=/home/user/scratch/cs581; R=/opt/work/iqd/reps; OUT=$H/iqtree/results_d/iqtree.jsonl
PY=/opt/mm/root/envs/pasta183/bin/python
for rep in "$@"; do
  echo "== $rep vote $(date)"
  python3 $H/gcmvote/code/run.py $R/$rep magus es4 hard-bb
  cp $R/$rep/results.jsonl $H/iqtree/results_d/$rep.vote.jsonl
  echo "== $rep iqtree $(date)"
  for m in true magus es4 hard-bb; do
    [ $m = true ] && a=$R/$rep/true.fasta || a=$R/$rep/vote/$m/out.fasta
    echo "$rep $m $a $R/$rep/true_tree.nwk $R/$rep/iq $OUT"
  done | xargs -P 2 -L 1 $PY $H/iqtree/code/iq_one.py
  echo "== $rep done $(date)"
done
