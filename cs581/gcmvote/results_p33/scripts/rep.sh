#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# rep.sh NAME : alignment lane + es4 + vote variants + 6 FastTree trees (4 at a time)
set -u
N=$1; lvl=${N%_R*}; r=${N#*_R}
S=/home/user/scratch/cs581
G=/opt/work/gcmtrees/reps/$N; V=/opt/work/gcmvote/reps/$N; T=/opt/work/gcmtrees/trees/${N}_d0
echo "[$(date +%T)] start $N"
bash $S/gcmtrees/code/run.sh $N
python3 $S/gcmgen/code/gg.py run $G 'linsi#es4'
echo "[$(date +%T)] gg done"
mkdir -p $V
for f in inputs true.fasta; do [ -e $V/$f ] || ln -s $G/$f $V/$f; done
python3 $S/gcmvote/code/run.py $V magus hard hard-bb soft4
echo "[$(date +%T)] vote done"
mkdir -p $T
ln -sf $G/true.fasta $T/true.fasta
ln -sf $G/variants/linsi/out.fasta $T/magus.fasta
ln -sf $G/variants/linsi_es_es4/out.fasta $T/es4.fasta
ln -sf $G/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta $T/recipe.fasta
ln -sf $V/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $V/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
ls -l $T
PY=/opt/mm/root/envs/pasta183/bin/python
mkdir -p /opt/work/tlog
pids=()
for m in true magus es4 recipe vote_hard vote_hard-bb; do
  while [ $(jobs -rp | wc -l) -ge 4 ]; do sleep 10; done
  (cd /tmp && nice -n 10 $PY $S/protbench/code/trees.py $T /opt/data/sim/$lvl/R$r/tree.nwk $G/trees.jsonl $m > /opt/work/tlog/${N}_$m.log 2>&1) &
  sleep 2
done
wait
echo "[$(date +%T)] REP DONE $N"
