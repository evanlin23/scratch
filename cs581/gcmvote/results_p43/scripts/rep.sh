#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# rep.sh NAME : one SIMHIGH replicate end to end (helper p43). Restartable: each stage skips done work.
#   gcmtrees run.sh (data + MAGUS draw + gcmtrees variants), gg.py linsi#es4, vote.py variants in a separate
#   rep dir, FastTree (4 independent trees.py processes at a time), bank tarball.
set -u
name=$1; lvl=${name%_R*}; r=${name#*_R}
S=/home/user/scratch/cs581
G=/opt/work/gcmtrees/reps/$name; V=/opt/work/gcmvote/reps/$name; T=/opt/work/gcmvote/trees/${name}_d0
echo "[$(date +%T)] start $name"
bash $S/gcmtrees/code/run.sh $name
python3 $S/gcmgen/code/gg.py run $G 'linsi#es4'
echo "[$(date +%T)] gcmtrees done"
mkdir -p $V $T
for f in inputs true.fasta unaligned.fasta; do [ -e $V/$f ] || ln -s $G/$f $V/$f; done
python3 $S/gcmvote/code/run.py $V magus hard hard-bb soft4
echo "[$(date +%T)] vote done"
ln -sf $G/true.fasta $T/true.fasta
ln -sf $G/variants/linsi/out.fasta $T/magus.fasta
ln -sf $G/variants/linsi_es_es4/out.fasta $T/es4.fasta
ln -sf $G/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta $T/recipe.fasta
ln -sf $V/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $V/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
for m in true magus es4 recipe vote_hard vote_hard-bb; do [ -s $T/$m.fasta ] || echo "MISSING $m"; done
printf '%s\n' true magus es4 recipe vote_hard vote_hard-bb | (cd /tmp && xargs -P 4 -I{} nice -n 10 \
  /opt/mm/root/envs/pasta183/bin/python $S/protbench/code/trees.py $T /opt/data/sim/$lvl/R$r/tree.nwk $V/trees.jsonl {})
echo "[$(date +%T)] trees done"
echo REP_DONE $name
