#!/bin/bash
# Restartable alignment lane: simulate, one MAGUS draw (paper flags), then the paired merge-only variants.
#   bash run.sh [DATASET ...]   (default: SIMHIGH_R1..R8 then SIMMOD_R1..R8)
set -u
G=/home/user/scratch/cs581/gcmtrees
C=/home/user/scratch/cs581/code
W=/opt/work/gcmtrees
V='linsi wsoft0.03:linsi&fftns2#es4 linsi#es3 linsi&fftns2-op3'
DS="$*"; [ -z "$DS" ] && DS="$(for l in SIMHIGH SIMMOD; do for r in 1 2 3 4 5 6 7 8; do echo ${l}_R$r; done; done)"
mkdir -p $W/reps
cd $C
for name in $DS; do
  lvl=${name%_R*}; r=${name#*_R}; draw=${DRAW:-0}
  bash $G/code/gen.sh $lvl $r >/dev/null || { echo "gen failed $name"; continue; }
  echo "$name /opt/data/sim/$lvl/R$r/sim.fa 25" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W --draws $draw --tools '' --e2e '' --union '' --threads 4 || continue
  rep=$W/reps/$name; [ "$draw" = 0 ] || rep=$W/reps/${name}_d$draw
  python3 /home/user/scratch/cs581/protcons/code/pc.py rep $W/${name}_d$draw $rep || continue
  python3 $G/code/../../gcmgen/code/gg.py run $rep $V
  bash $G/code/collect.sh
done
echo ALN_DONE
