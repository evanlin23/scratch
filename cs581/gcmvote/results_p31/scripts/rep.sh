#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# rep.sh NAME : alignment lane + vote variants + 6 FastTree trees (4 parallel) for one replicate
set -u
n=$1; lvl=${n%_R*}; r=${n#*_R}
S=/home/user/scratch/cs581
R=/opt/work/gcmtrees/reps/$n; V=/opt/work/gcmvote/reps/$n
echo "== aln $(date)"
bash $S/gcmtrees/code/run.sh $n
python3 $S/gcmgen/code/gg.py run $R 'linsi#es4'
echo "== vote $(date)"
mkdir -p $V; ln -sfn $R/inputs $V/inputs; ln -sf $R/true.fasta $V/true.fasta
python3 $S/gcmvote/code/run.py $V magus hard hard-bb soft4
echo "== trees $(date)"
T=/opt/work/gcmtrees/trees/${n}_d0; mkdir -p $T
ln -sf $R/true.fasta $T/true.fasta
ln -sf $R/variants/linsi/out.fasta $T/magus.fasta
ln -sf $R/variants/linsi_es_es4/out.fasta $T/es4.fasta
ln -sf $R/variants/wsoft0.03_c_linsi_i_fftns2_es_es4/out.fasta $T/recipe.fasta
ln -sf $V/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $V/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
cd /tmp
printf '%s\n' true magus es4 recipe vote_hard vote_hard-bb | xargs -P 4 -I{} nice -n 10 /opt/mm/root/envs/pasta183/bin/python $S/protbench/code/trees.py $T /opt/data/sim/$lvl/R$r/tree.nwk $R/trees.jsonl {}
echo "== REP_DONE $n $(date)"
