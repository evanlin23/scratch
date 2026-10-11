#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# rep.sh NAME : after gcmtrees run.sh + gg.py 'linsi#es4': vote variants in a separate rep dir, then
# 6 FastTree runs (4 at a time, independent processes). Restartable.
set -u
n=$1; G=/opt/work/gcmtrees/reps/$n; V=/opt/work/gcmvote/reps/$n; H=$(dirname $(readlink -f $0))
mkdir -p $V; ln -sfn $G/inputs $V/inputs; ln -sf $G/true.fasta $V/true.fasta
python3 /home/user/scratch/cs581/gcmvote/code/run.py $V magus hard hard-bb soft4 || exit 1
printf '%s\n' "true $G/true.fasta" "magus $G/variants/linsi/out.fasta" "es4 $G/variants/linsi_es_es4/out.fasta" \
  "recipe $G/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta" "vote_hard $V/vote/hard/out.fasta" \
  "vote_hard-bb $V/vote/hard-bb/out.fasta" | xargs -P 4 -L 1 bash $H/tree1.sh $n
echo "REP_DONE $n"
