#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# rep.sh K : SIMHIGH_RK end to end (alignment lane, vote variants, 6 FastTree trees 4 at a time). Restartable.
set -u
K=$1; N=SIMHIGH_R$K; C=/home/user/scratch/cs581
R=/opt/work/gcmtrees/reps/$N; V=/opt/work/gcmvote/reps/$N; T=/opt/work/gcmtrees/trees/${N}_d0
echo "[$(date +%T)] start $N"
bash $C/gcmtrees/code/run.sh $N
python3 $C/gcmgen/code/gg.py run $R 'linsi#es4'
echo "[$(date +%T)] gcmtrees lane done"
mkdir -p $V
for f in inputs true.fasta; do [ -e $V/$f ] || ln -s $R/$f $V/$f; done
python3 $C/gcmvote/code/run.py $V magus hard hard-bb soft4
echo "[$(date +%T)] vote done"
mkdir -p $T
ln -sf $R/true.fasta $T/true.fasta
ln -sf $R/variants/linsi/out.fasta $T/magus.fasta
ln -sf $R/variants/linsi_es_es4/out.fasta $T/es4.fasta
ln -sf $R/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta $T/recipe.fasta
ln -sf $V/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $V/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
mkdir -p $T/rows
printf '%s\n' true magus es4 recipe vote_hard vote_hard-bb | xargs -P 4 -I{} sh -c \
  "cd /tmp && nice -n 10 /opt/mm/root/envs/pasta183/bin/python $C/protbench/code/trees.py $T /opt/data/sim/SIMHIGH/R$K/tree.nwk $T/rows/{}.jsonl {} > $T/rows/{}.log 2>&1"
cat $T/rows/*.jsonl > $R/trees.jsonl
echo "[$(date +%T)] REP_DONE $N"
