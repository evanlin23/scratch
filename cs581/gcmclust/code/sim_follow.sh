#!/bin/bash
# Prep + build each AliSim rep once its fresh MAGUS draw is recorded (restartable).
C=/home/user/scratch/cs581/gcmclust/code; W=/opt/work/gcmclust
for n in SIMMOD_R1 SIMHIGH_R1 SIMMOD_R2 SIMHIGH_R2; do
  until grep -q "\"$n\"" $W/fresh_magus.jsonl 2>/dev/null; do sleep 60; done
  k=${n%_R*}; r=${n#*_R}
  python3 $C/gc.py prep $n $W/reps/$n /opt/data/sim/$k/R$r/sim.fa $W/fresh/${n}_d0 && python3 $C/gc.py build $W/reps/$n
done
echo SIMFOLLOW_DONE
