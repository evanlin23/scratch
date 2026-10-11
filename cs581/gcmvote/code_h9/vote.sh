#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Pre-registered vote variants on every h9 rep (after sens.sh). Restartable.
W=/opt/work/h9; H=/home/user/scratch/cs581/gcmvote/code_h9
for name in 1000M2_R0 1000L1_R0 16S.M_R0; do
  for B in 5 10 20; do
    [ -d $W/reps/${name}_B$B/sets/s0 ] && python3 $H/vote_run.py $W/reps/${name}_B$B magus hard hard-bb soft soft-bb soft2 soft4
  done
done
echo VOTE_DONE
