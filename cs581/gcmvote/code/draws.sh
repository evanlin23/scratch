#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Restartable: AliSim data (gcmtrees gen.sh = protbench seeds), one fresh MAGUS draw (paper flags), replicate dir.
#   bash draws.sh NAME [NAME ...]   e.g. SIMHIGH_R1
C=/home/user/scratch/cs581/code; W=/opt/work/gcmvote
cd $C
for name in "$@"; do
  rep=$W/reps/$name; [ -d $rep/sets/s0 ] && continue
  lvl=${name%_R*}; r=${name#*_R}
  bash /home/user/scratch/cs581/gcmtrees/code/gen.sh $lvl $r >/dev/null || { echo "gen failed $name"; continue; }
  echo "$name /opt/data/sim/$lvl/R$r/sim.fa 25" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/draws --draws 0 --tools '' --e2e '' --union '' --threads 4 || continue
  python3 /home/user/scratch/cs581/protcons/code/pc.py rep $W/draws/${name}_d0 $rep
done
echo DRAWS_DONE
