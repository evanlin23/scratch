#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# pipe.sh NAME (from helper p25): gcmtrees lane + gg es4 + vote variants (separate rep dir) + 6 FastTree trees (4 at a time)
set -u
N=$1; lvl=${N%_R*}; r=${N#*_R}
S=/home/user/scratch/cs581; R=/opt/work/gcmtrees/reps/$N; V=/opt/work/gcmvote/reps/$N
T=/opt/work/gcmtrees/trees/${N}_d0
echo "== start $N $(date)"
bash $S/gcmtrees/code/run.sh $N
python3 $S/gcmgen/code/gg.py run $R 'linsi#es4'
echo "== aln done $N $(date)"
mkdir -p $V; ln -sfn $R/inputs $V/inputs; ln -sf $R/true.fasta $V/true.fasta
python3 $S/gcmvote/code/run.py $V magus hard hard-bb soft4
echo "== vote done $N $(date)"
mkdir -p $T
ln -sf $R/true.fasta $T/true.fasta
ln -sf $R/variants/linsi/out.fasta $T/magus.fasta
ln -sf $R/variants/linsi_es_es4/out.fasta $T/es4.fasta
ln -sf $R/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta $T/recipe.fasta
ln -sf $V/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $V/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
ls -l $T
printf '%s\n' true magus es4 recipe vote_hard vote_hard-bb | xargs -P 4 -I{} sh -c \
  "cd /tmp && nice -n 10 /opt/mm/root/envs/pasta183/bin/python $S/protbench/code/trees.py $T /opt/data/sim/$lvl/R$r/tree.nwk $R/trees.jsonl {} > /opt/work/tree_${N}_{}.log 2>&1"
echo "== PIPE DONE $N $(date)"
