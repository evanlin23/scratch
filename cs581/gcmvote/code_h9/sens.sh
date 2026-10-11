#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Backbone-count sensitivity (helper h9): one MAGUS draw with 20 backbones per dataset, then merge-only
# reps on the same subsets with B = 5, 10, 20 (backbones 1..B), MAGUS's merge and edge-support filters.
# Restartable: bbtool_bench resumes from state.json, gg.py skips variants already in results.jsonl.
#   bash sens.sh [DATASET ...]
W=/opt/work/h9; R=/home/user/scratch; C=$R/cs581/code; G=$R/cs581/gcmgen/code
mkdir -p $W
cd $C
[ $# -eq 0 ] && set -- 1000M2_R0 1000L1_R0 16S.M_R0
for name in "$@"; do
  grep -h "^$name " fanout/jobs_*.txt | head -1 | awk '{print $1, $2, $3}' > $W/job_$name.txt
  GCMX_NUM_BACKBONES=20 python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/fresh --draws 0 \
    --tools '' --e2e '' --union '' --threads 4 || continue
  python3 $R/cs581/protcons/code/pc.py rep $W/fresh/${name}_d0 $W/reps/${name}_B20 || continue
  for B in 5 10; do
    d=$W/reps/${name}_B$B
    if [ ! -d $d/sets/s0 ]; then
      rm -rf $d; mkdir -p $d/inputs/backbones
      cp -r $W/reps/${name}_B20/inputs/subalignments $d/inputs/
      for i in $(seq 1 $B); do cp $W/reps/${name}_B20/inputs/backbones/backbone_${i}_mafft.txt $d/inputs/backbones/; done
      cp $W/reps/${name}_B20/{unaligned.fasta,magus.json} $d/
      python3 -c "import sys; sys.path.insert(0,'$R/cs581/bbevidence/code'); import bbe; bbe.prep('$d', '$W/reps/${name}_B20/true.fasta')"
    fi
  done
  for B in 5 10 20; do
    v="linsi"
    for f in 0.2 0.3 0.4 0.5; do v="$v linsi#es$(python3 -c "import math; print(math.ceil($f*$B-1e-9))")"; done
    python3 $G/gg.py run $W/reps/${name}_B$B $(echo $v | tr ' ' '\n' | awk '!s[$0]++')
  done
done
echo SENS_DONE
