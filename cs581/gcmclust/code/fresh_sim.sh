#!/bin/bash
# AliSim protein replicates (as cs581/protbench, same seeds) + one fresh MAGUS draw each (paper flags). Restartable.
C=/home/user/scratch/cs581/gcmclust/code; W=/opt/work/gcmclust
for j in "MOD 1" "HIGH 1" "MOD 2" "HIGH 2"; do set -- $j; bash $C/sim.sh $1 $2; done
cd /home/user/scratch/cs581/code
for j in SIMMOD_R1 SIMHIGH_R1 SIMMOD_R2 SIMHIGH_R2; do
  k=${j%_R*}; r=${j#*_R}
  echo "$j /opt/data/sim/$k/R$r/sim.fa 25" > $W/job_$j.txt
  python3 -m gcmx.bbtool_bench $W/job_$j.txt $W/fresh_magus.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads 4
done
echo FRESH_DONE
