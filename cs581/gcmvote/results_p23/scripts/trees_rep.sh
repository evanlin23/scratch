#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# trees_rep.sh SIMHIGH_Rk : FastTree -lg -gamma via protbench trees.py, 4 independent processes at a time
set -u
N=$1; R=/opt/work/gcmtrees/reps/$N; V=/opt/work/gcmvote/reps/$N; T=/opt/work/gcmtrees/trees/${N}_d0
r=${N#*_R}; mkdir -p $T
ln -sf $R/true.fasta $T/true.fasta
ln -sf $R/variants/linsi/out.fasta $T/magus.fasta
ln -sf $R/variants/linsi_es_es4/out.fasta $T/es4.fasta
ln -sf $R/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta $T/recipe.fasta
ln -sf $V/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $V/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
export T R r
printf '%s\n' true magus es4 recipe vote_hard vote_hard-bb | xargs -P 4 -I{} bash -c \
  'cd /tmp && nice -n 10 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T /opt/data/sim/SIMHIGH/R$r/tree.nwk $R/trees.jsonl {} > /opt/work/tree_${r}_{}.log 2>&1'
echo "TREES_DONE $N"
