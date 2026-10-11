#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# rep.sh K : SIMHIGH_RK alignment lane + gg es4 + vote variants + 6 FastTree trees (4 at a time)
set -u
K=$1; N=SIMHIGH_R$K; S=/home/user/scratch/cs581
R=/opt/work/gcmtrees/reps/$N; V=/opt/work/gcmvote/reps/$N; T=/opt/work/gcmtrees/trees/${N}_d0
echo "[$(date +%T)] start $N"
bash $S/gcmtrees/code/run.sh $N
echo "[$(date +%T)] run.sh done"
python3 $S/gcmgen/code/gg.py run $R 'linsi#es4'
echo "[$(date +%T)] gg es4 done"
mkdir -p $V; ln -sfn $R/inputs $V/inputs; ln -sf $R/true.fasta $V/true.fasta
python3 $S/gcmvote/code/run.py $V magus hard hard-bb soft4
echo "[$(date +%T)] vote done"
mkdir -p $T
ln -sf $R/true.fasta $T/true.fasta
ln -sf $R/variants/linsi/out.fasta $T/magus.fasta
ln -sf $R/variants/linsi_es_es4/out.fasta $T/es4.fasta
ln -sf $R/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta $T/recipe.fasta
ln -sf $V/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $V/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
ls -l $T
cd /tmp
printf '%s\n' true magus es4 recipe vote_hard vote_hard-bb | xargs -P 4 -I{} nice -n 10 /opt/mm/root/envs/pasta183/bin/python $S/protbench/code/trees.py $T /opt/data/sim/SIMHIGH/R$K/tree.nwk $R/trees.jsonl {}
echo "[$(date +%T)] trees done"
echo REP_DONE $N
