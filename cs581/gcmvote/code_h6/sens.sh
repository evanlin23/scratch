#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Restartable: one MAGUS draw with 20 backbones per dataset, then merge-only variants for B = 5, 10, 20
# (the first B backbones: bbe variant linsi~B, vote.py --bb REP/bbB) -> REP/{results,vote_results}.jsonl;
# collect.py -> results_h6/sens.jsonl. Post-MAGUS merges run in parallel lanes (merge_wall is under contention).
#   bash sens.sh DATASET [...]
set -u
H=/home/user/scratch/cs581/gcmvote
C=/home/user/scratch/cs581/code
W=/opt/work/h6
mkdir -p $W/reps
cd $C; export PYTHONPATH=$C
V='linsi~5 linsi~5#es1 linsi~5#es2 linsi~5#es3 linsi~10 linsi~10#es2 linsi~10#es3 linsi~10#es4 linsi~10#es5 linsi linsi#es4 linsi#es6 linsi#es8 linsi#es10'
for name in "$@"; do
  case $name in
    BBA*) src=$C/../data/balibase_clean/RV100_$name.fasta;;
    *) lvl=${name%_R*}; r=${name#*_R}; bash $H/../gcmtrees/code/gen.sh $lvl $r >/dev/null || { echo "gen failed $name"; continue; }
       src=/opt/data/sim/$lvl/R$r/sim.fa;;
  esac
  echo "$name $src 25" > $W/job_$name.txt
  python3 $H/code_h6/magus20.py $W/job_$name.txt $W/magus.jsonl $W --draws 0 --tools '' --e2e '' --union '' --threads 4 || continue
  rep=$W/reps/$name
  python3 $H/../protcons/code/pc.py rep $W/${name}_d0 $rep || continue
  rm -f $rep/variants/*.lock  # stale locks from killed lanes (gg.py skips locked variants)
  python3 $H/../gcmgen/code/gg.py run $rep $V > $rep/lane_gg1.log 2>&1 &
  sleep 2; python3 $H/../gcmgen/code/gg.py run $rep $(echo $V | tr ' ' '\n' | tac) > $rep/lane_gg2.log 2>&1 &
  for B in 5 10 20; do python3 $H/code_h6/vrun.py $rep $B hard soft soft2 soft4 hard-bb soft-bb > $rep/lane_v$B.log 2>&1 & done
  wait
  python3 $H/code_h6/collect.py
done
echo SENS_DONE
