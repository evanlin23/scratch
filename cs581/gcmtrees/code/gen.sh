#!/bin/bash
# AliSim protein replicates exactly as cs581/protbench (REPORT.md "Data"): tree seed 100r+7, sequence seed 100r+13.
#   bash gen.sh LEVEL R   (LEVEL = SIMMOD | SIMHIGH | SIMVHIGH)  -> /opt/data/sim/LEVEL/R<r>/{tree.nwk,sim.fa}
set -euo pipefail
lvl=$1; r=$2
case $lvl in SIMMOD) mean=0.06;; SIMHIGH) mean=0.10;; SIMVHIGH) mean=0.14;; esac
IQ=/opt/mm/root/envs/bio/bin/iqtree3
D=/opt/data/sim/$lvl/R$r; mkdir -p $D; cd $D
[ -s sim.fa ] && exit 0
$IQ -r 1000 tree.nwk -rlen 0.001 $mean 0.8 -seed $((100*r+7)) -redo > gen_tree.log 2>&1
$IQ --alisim sim -t tree.nwk -m LG+G4 --length 300 --indel 0.05,0.05 \
    --indel-size POW{1.7/40},POW{1.7/40} -seed $((100*r+13)) -af fasta -redo > gen_seq.log 2>&1
ls
