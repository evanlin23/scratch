#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# FastTree -lg -gamma (protbench trees.py) on the true alignment and vote2 outputs of a SIMHIGH bank rep;
# nRF vs the bank's true tree -> REP/trees2.jsonl. Restartable.   bash trees2.sh REPDIR VARIANT [...]
rep=$1; shift; name=$(basename $rep)
T=/opt/work/trees/${name}_d0; mkdir -p $T
ln -sf $rep/true.fasta $T/true.fasta
ms="true"
for v in "$@"; do src=$rep/vote2/$v/out.fasta; [ -s $src ] || { echo "missing $src"; continue; }; ln -sf $src $T/$v.fasta; ms="$ms $v"; done
nice -n 5 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T $rep/true_tree.nwk $rep/trees2.jsonl $ms
