#!/bin/bash
# TreeMerge on the simulated replicates of the three pilot conditions (2 in parallel).
cd "$(dirname "$0")"
for c in "yule_n200_m50_i0.005 exact" "yule_n200_m50_i0.005 iqtree" "cat_n200_m50_i0.002 exact"; do
  read -r S ST <<< "$c"
  seq 1 20 | xargs -P 2 -I{} python3 sim_treemerge.py /opt/gtmdata/simq/$S/{} $ST
done
echo TMDONE
