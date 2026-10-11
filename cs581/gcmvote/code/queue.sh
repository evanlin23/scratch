#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Restartable lane: lines "REP B VARIANT [VARIANT ...]"; reps under /opt/work/gcmvote/reps; skips missing reps.
G=/home/user/scratch/cs581/gcmvote/code
while read -r rep B rest; do
  [ -z "$rep" ] && continue; [[ $rep == \#* ]] && continue
  [ -d /opt/work/gcmvote/reps/$rep/sets/s0 ] || { echo "skip $rep (no rep yet)"; continue; }
  python3 $G/run.py /opt/work/gcmvote/reps/$rep --B $B $rest
done < "$1"
echo QUEUE_DONE "$1"
