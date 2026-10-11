#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# tree1.sh NAME METHOD SRC : link SRC as METHOD.fasta, run FastTree -lg -gamma via protbench trees.py
n=$1; T=/opt/work/gcmtrees/trees/${n}_d0; R=/opt/work/gcmtrees/reps/$n; r=${n#*_R}
mkdir -p $T; ln -sf $3 $T/$2.fasta
cd /tmp && nice -n 10 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py \
  $T /opt/data/sim/SIMHIGH/R$r/tree.nwk $R/trees.jsonl $2 > /opt/work/tree_${n}_$2.log 2>&1
