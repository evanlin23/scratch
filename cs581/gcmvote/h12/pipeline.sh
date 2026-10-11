#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# After the MAGUS draw: gg.py baselines, then every pre-declared vote variant (B = 10), each under memwatch.
# Vote rows go to a second rep dir (symlinked inputs) because run.py's results.jsonl rows carry "B".
W=/opt/work/h12; H=/home/user/scratch/cs581/gcmvote/h12; REP=$W/reps/16S.T_R0; V=$W/reps/16S.T_R0v
VCODE=$W/vcode/cs581/gcmvote/code
until grep -q MAGUS_DONE $W/magus.out; do grep -q MAGUS_FAILED $W/magus.out && exit 1; sleep 30; done
bash $H/variants.sh gg linsi 'linsi#es3' 'linsi#es4' 'linsi#es5'
mkdir -p $V
for f in inputs true.fasta unaligned.fasta magus.json sets; do [ -e $V/$f ] || ln -s $REP/$f $V/$f; done
cd /home/user/scratch/cs581/code
mkdir -p $W/mem
for v in magus es4 hard hard-bb soft soft-bb soft2 soft4 frac0.2 frac0.3 frac0.4 frac0.5 hard+mask; do
  [ -f $V/results.jsonl ] && grep -q "\"variant\": \"$v\", \"B\": 10" $V/results.jsonl && continue
  python3 $H/memwatch.py $W/mem/vote_$v.json 14400 python3 $VCODE/run.py $V $v || echo "VOTE_FAILED $v"
done
echo PIPELINE_DONE
