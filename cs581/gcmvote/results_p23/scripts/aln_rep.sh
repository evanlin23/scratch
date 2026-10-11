#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# aln_rep.sh SIMHIGH_Rk : gcmtrees lane + gg.py es4 + vote variants in a separate rep dir
set -u
N=$1; S=/home/user/scratch/cs581
bash $S/gcmtrees/code/run.sh $N
R=/opt/work/gcmtrees/reps/$N
python3 $S/gcmgen/code/gg.py run $R 'linsi#es4'
V=/opt/work/gcmvote/reps/$N; mkdir -p $V
for f in inputs true.fasta unaligned.fasta; do [ -e $V/$f ] || ln -s $R/$f $V/$f; done
python3 $S/gcmvote/code/run.py $V magus hard hard-bb soft4
echo "REP_DONE $N"
