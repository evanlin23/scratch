#!/bin/bash
# AliSim protein replicates exactly as cs581/protbench (REPORT "Simulation"): Yule-Harding tree on 1,000 taxa,
# LG+G4, root length 300, indels 0.05/0.05 with POW{1.7/40} sizes; tree seed 100r+7, sequence seed 100r+13.
# Usage: sim.sh MOD|HIGH R
set -eu
IQ=/opt/mm/root/envs/bio/bin/iqtree3
kind=$1; r=$2
mean=$([ "$kind" = MOD ] && echo 0.06 || echo 0.10)
D=/opt/data/sim/SIM$kind/R$r
[ -f $D/sim.fa ] && exit 0
mkdir -p $D; cd $D
$IQ -r 1000 -rlen 0.001 $mean 0.8 tree.nwk -seed $((100*r+7)) -redo > tree.log 2>&1
$IQ --alisim sim -t tree.nwk -m LG+G4 --length 300 --indel 0.05,0.05 \
  --indel-size "POW{1.7/40},POW{1.7/40}" -af fasta -seed $((100*r+13)) -redo > alisim.log 2>&1
ls $D
