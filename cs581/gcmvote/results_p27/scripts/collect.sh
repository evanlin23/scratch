#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# rebuild results_p27 from the work dirs
O=/home/user/scratch/cs581/gcmvote/results_p27; mkdir -p $O/scripts
: > $O/aln.jsonl; : > $O/trees.jsonl
for n in SIMHIGH_R27 SIMHIGH_R28; do
  R=/opt/work/gcmtrees/reps/$n; V=/opt/work/gcmvote/reps/$n
  [ -f $R/results.jsonl ] && cat $R/results.jsonl >> $O/aln.jsonl
  [ -f $V/results.jsonl ] && cat $V/results.jsonl >> $O/aln.jsonl
  for m in true magus es4 recipe vote_hard vote_hard-bb; do [ -f $R/trees.$m.jsonl ] && cat $R/trees.$m.jsonl >> $O/trees.jsonl; done
done
cp /opt/work/gcmtrees/magus.jsonl $O/magus.jsonl
cp $(dirname $0)/pipe.sh $(dirname $0)/bank.sh $(dirname $0)/norm.py $0 $O/scripts/
