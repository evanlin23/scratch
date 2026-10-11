#!/bin/bash
# Regenerate AliSim SIMHIGH replicates exactly as claude/cs581-protbench (REPORT.md "Simulation"):
# Yule-Harding tree on 1000 taxa (iqtree3 -r 1000 -rlen 0.001 0.10 0.8, seed 100r+7), then
# LG+G4, root length 300, indels 0.05/0.05, POW{1.7/40} sizes (seed 100r+13).
set -euo pipefail
IQ=/opt/mm/root/envs/bio/bin/iqtree3
for r in "$@"; do
  d=/opt/data/sim/SIMHIGH/R$r; mkdir -p $d; cd $d
  [ -f sim.fa ] && continue
  $IQ -r 1000 -rlen 0.001 0.10 0.8 tree.nwk -seed $((100*r+7)) -redo -quiet > /dev/null 2>&1 || true  # -redo: iqtree3 touches the output before checking it
  $IQ --alisim sim -t tree.nwk -m LG+G4 --length 300 --indel 0.05,0.05 --indel-size "POW{1.7/40},POW{1.7/40}" \
      -seed $((100*r+13)) -af fasta -quiet > /dev/null
  ls
done
