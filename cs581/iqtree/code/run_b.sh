#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Per SIMHIGH replicate: rebuild magus / es4 / hard-bb merges with gcmvote's run.py (scored with FastSP),
# then IQ-TREE (iq_b.py) on true + the three merges, 2 runs at a time. Restartable.
#   bash run_b.sh SIMHIGH_R3 [SIMHIGH_R4 ...]      (reps unpacked under /opt/work/iqb/reps)
set -u
S=/home/user/scratch/cs581; W=/opt/work/iqb; OUT=$S/iqtree/results_b/iqtree.jsonl
PY=/opt/mm/root/envs/pasta183/bin/python
for name in "$@"; do
  rep=$W/reps/$name
  python3 $S/gcmvote/code/run.py $rep magus es4 hard-bb || { echo "MERGE FAIL $name"; continue; }
  cp $rep/results.jsonl $S/iqtree/results_b/$name.merge.jsonl
  printf '%s\n' "true $rep/true.fasta" "magus $rep/vote/magus/out.fasta" "es4 $rep/vote/es4/out.fasta" \
    "hard-bb $rep/vote/hard-bb/out.fasta" | \
    xargs -P 2 -L 1 bash -c '$0 '$S'/iqtree/code/iq_b.py '$name' $1 $2 '$rep'/true_tree.nwk '$W'/iq/'$name' '$OUT' \
      || echo "IQ FAIL '$name' $1"' $PY
  echo "REP DONE $name"
done
echo "ALL DONE"
