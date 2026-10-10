#!/bin/bash
# IQ-TREE lane over SIMHIGH replicates whose FastTree trees exist (true, magus, recipe first; then es3, hard).
#   bash iqtrees.sh [PARALLEL] [METHODS...]
W=/opt/work/gcmtrees; G=/home/user/scratch/cs581/gcmtrees
P=${1:-2}; shift; M="${*:-true magus recipe}"
one() {
  T=$1; n=$(basename $T); n=${n%_d0}; b=${n%_d[0-9]*}; lvl=${b%_R*}; r=${b#*_R}
  nice -n 5 /opt/mm/root/envs/pasta183/bin/python $G/code/iqtrees.py $T /opt/data/sim/$lvl/R$r/tree.nwk $W/reps/$n/iqtrees.jsonl $M >/dev/null
}
export -f one; export W G M
ls -d $W/trees/SIMHIGH_* 2>/dev/null | xargs -r -P $P -I{} bash -c 'one {}'
bash $G/code/collect.sh
echo IQ_DONE
