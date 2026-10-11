#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Calibrate every replicate whose magus vote run exists (restartable); loops until STOP_CALIB exists.
W=/opt/work/gcmvote
while [ ! -f $W/STOP_CALIB ]; do
  for d in $W/reps/*/; do
    [ -f $d/vote/magus/edges.npz ] && [ -f $d/vote/magus/run.json ] && [ ! -f $d/calib.json ] && \
      nice -n 5 python3 /home/user/scratch/cs581/gcmvote/code/calib.py $d > /dev/null 2>$d/calib.err && echo "calib $d"
  done
  sleep 60
done
