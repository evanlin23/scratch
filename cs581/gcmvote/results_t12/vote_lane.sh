#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# vote variants for SIMHIGH_R12 on a separate rep dir (results.jsonl formats differ from gg.py's)
R=/opt/work/gcmtrees/reps/SIMHIGH_R12; V=/opt/work/gcmvote/reps/SIMHIGH_R12
until [ -s $R/true.fasta ] && [ -d $R/sets/s0 ]; do sleep 20; done
mkdir -p $V
[ -e $V/inputs ] || ln -s $R/inputs $V/inputs
[ -e $V/true.fasta ] || ln -s $R/true.fasta $V/true.fasta
[ -e $V/sets ] || ln -s $R/sets $V/sets
cd /home/user/scratch
printf '%s\n' magus es4 hard soft soft2 soft4 hard-bb soft-bb hard+mask | \
  xargs -P ${P:-3} -I{} python3 cs581/gcmvote/code/run.py $V {}
echo VOTE_DONE
