#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# tree1.sh METHOD SRC  : link SRC as METHOD.fasta, run FastTree via trees.py
T=/opt/work/gcmtrees/trees/SIMHIGH_R20_d0; R=/opt/work/gcmtrees/reps/SIMHIGH_R20
ln -sf $2 $T/$1.fasta
cd /tmp && nice -n 10 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T /opt/data/sim/SIMHIGH/R20/tree.nwk $R/trees.jsonl $1 > /opt/work/tree_$1.log 2>&1
