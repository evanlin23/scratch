#!/bin/bash
# Post-hoc baselines requested after the pre-registration was frozen (REPORT.md §post-hoc):
# global edge-support threshold k in {2,3,5} (run_edgesup.py). Optional: VARS env overrides the list
# (e.g. VARS='linsi@n20' for 20 L-INS-i backbones). Restartable.
P=/home/user/scratch/cs581/protcons
W=/opt/work/protcons
V=${VARS:-'linsi#es2 linsi#es3 linsi#es5'}
for R in $W/reps/*/; do
  n=$(basename $R)
  [ -f $R/gate.json ] || continue
  [ -f $R/unaligned.fasta ] || { [ -f $W/${n}_d0/unaligned.fasta ] && cp $W/${n}_d0/unaligned.fasta $R/; }
  python3 $P/code/pc.py run $R $V
  bash $P/code/collect.sh
done
