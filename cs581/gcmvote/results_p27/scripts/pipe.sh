#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# pipe.sh SIMHIGH_Rk : gcmtrees lane + linsi#es4 + vote variants + 6 FastTree trees (4 at a time)
set -u
n=$1; lvl=${n%_R*}; r=${n#*_R}
C=/home/user/scratch/cs581; W=/opt/work/gcmtrees; R=$W/reps/$n; V=/opt/work/gcmvote/reps/$n
echo "[$(date +%T)] start $n"
bash $C/gcmtrees/code/run.sh $n
python3 $C/gcmgen/code/gg.py run $R 'linsi#es4'
echo "[$(date +%T)] gg done $n"
mkdir -p $V; ln -sfn $R/inputs $V/inputs; ln -sf $R/true.fasta $V/true.fasta
python3 $C/gcmvote/code/run.py $V magus hard hard-bb soft4
echo "[$(date +%T)] vote done $n"
T=$W/trees/${n}_d0; mkdir -p $T
ln -sf $R/true.fasta $T/true.fasta
ln -sf $R/variants/linsi/out.fasta $T/magus.fasta
ln -sf $R/variants/linsi_es_es4/out.fasta $T/es4.fasta
ln -sf $R/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta $T/recipe.fasta
ln -sf $V/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $V/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
ls -lL $T
cd /tmp
printf '%s\n' true magus es4 recipe vote_hard vote_hard-bb | xargs -P 4 -I{} sh -c \
  "nice -n 10 /opt/mm/root/envs/pasta183/bin/python $C/protbench/code/trees.py $T /opt/data/sim/$lvl/R$r/tree.nwk $R/trees.{}.jsonl {} > /opt/work/tree_${n}_{}.log 2>&1"
cat $R/trees.true.jsonl $R/trees.magus.jsonl $R/trees.es4.jsonl $R/trees.recipe.jsonl $R/trees.vote_hard.jsonl $R/trees.vote_hard-bb.jsonl > $R/trees_p27.jsonl
echo "[$(date +%T)] PIPE_DONE $n"
