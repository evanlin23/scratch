#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# trees.sh NAME : FastTree (trees.py) on true, magus, es4, recipe, vote_hard, vote_hard-bb; 4 independent processes
n=$1; lvl=${n%_R*}; r=${n#*_R}
G=/opt/work/gcmtrees/reps/$n; V=/opt/work/gcmvote/reps/$n; T=/opt/work/gcmvote/trees/${n}_d0; mkdir -p $T
ln -sf /opt/work/gcmtrees/${n}_d0/true.fasta $T/true.fasta
ln -sf /opt/work/gcmtrees/${n}_d0/magus.fasta $T/magus.fasta
ln -sf $G/variants/$ES4/out.fasta $T/es4.fasta
ln -sf $G/variants/$RECIPE/out.fasta $T/recipe.fasta
ln -sf $V/vote/hard/out.fasta $T/vote_hard.fasta
ln -sf $V/vote/hard-bb/out.fasta $T/vote_hard-bb.fasta
for m in true magus es4 recipe vote_hard vote_hard-bb; do [ -s $T/$m.fasta ] || { echo MISSING $m; exit 1; }; done
cd /tmp
printf '%s\n' true magus es4 recipe vote_hard vote_hard-bb | xargs -P 4 -I{} sh -c \
  "nice -n 5 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T /opt/data/sim/$lvl/R$r/tree.nwk $V/trees.jsonl {} > $T/{}.trees.log 2>&1"
cat $V/trees.jsonl; echo "TREES_DONE $n"
