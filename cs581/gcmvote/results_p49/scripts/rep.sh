#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# rep.sh SIMHIGH_Rk : data + MAGUS draw + gcmtrees variants, gg.py linsi#es4, vote variants (separate rep dir),
# then FastTree (4 independent processes at a time) on true/magus/es4/recipe/vote_hard/vote_hard-bb. Restartable.
set -u
N=$1; lvl=${N%_R*}; r=${N#*_R}
S=/home/user/scratch/cs581
R=/opt/work/gcmtrees/reps/$N; VR=/opt/work/gcmvote/reps/$N; T=/opt/work/gcmtrees/trees/${N}_d0
echo "== $N start $(date)"
[ -s $R/variants/linsi_i_fftns2-op3/out.fasta ] || bash $S/gcmtrees/code/run.sh $N
python3 $S/gcmgen/code/gg.py run $R 'linsi#es4'
echo "== $N gg done $(date)"
mkdir -p $VR
for f in inputs true.fasta unaligned.fasta; do [ -e $R/$f ] && ln -sfn $R/$f $VR/$f; done
python3 $S/gcmvote/code/run.py $VR magus hard hard-bb soft4
echo "== $N vote done $(date)"
mkdir -p $T
ln -sf $R/true.fasta $T/true.fasta
ln -sf $R/variants/linsi/out.fasta $T/magus.fasta
ln -sf $R/variants/linsi_es_es4/out.fasta $T/es4.fasta
ln -sf $R/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta $T/recipe.fasta
ln -sf $VR/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $VR/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
for m in true magus es4 recipe vote_hard vote_hard-bb; do echo $m; done | xargs -P 4 -I{} bash -c \
  "cd /tmp && nice -n 10 /opt/mm/root/envs/pasta183/bin/python $S/protbench/code/trees.py $T /opt/data/sim/$lvl/R$r/tree.nwk $R/trees.jsonl {} > /opt/work/tree_${N}_{}.log 2>&1"
echo "== $N trees done $(date)"
echo "REP_DONE $N"
