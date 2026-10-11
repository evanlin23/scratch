#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Restartable backbone generation lane: per dataset s100 (10), s400 (10), s200 #21..40 (20).
# Seeds: SEEDBASE = 100000*dataset_index + 100*size (+ backbone index i).
C=/home/user/scratch/cs581/bbsize/code; W=/opt/work/bbsize
DS="BBA0101_b20 BBA0067_b20 SIMHIGH_R1 SIMMOD_R1 1000M2_R0_B20 16S.M_R0_B20"
i=0
for d in $DS; do
  i=$((i+1))
  python3 $C/bbgen.py $W/bank/$d $W/bb/$d/s100 100 1 10 $((100000*i+10000)) || echo "FAIL $d s100"
  python3 $C/bbgen.py $W/bank/$d $W/bb/$d/s400 400 1 10 $((100000*i+40000)) || echo "FAIL $d s400"
  python3 $C/bbgen.py $W/bank/$d $W/bb/$d/s200 200 21 40 $((100000*i+20000)) || echo "FAIL $d s200"
  touch $W/bb/$d/GEN_DONE
done
echo GEN_ALL_DONE
