#!/bin/bash
# Overnight queue (sequential, so wall-clock times are not confounded by other runs).
C=/home/user/scratch/cs581/epangdown/code
for D in /opt/work/nt78 /opt/work/16S; do
  $C/run_bscampp.sh $D stock 2000 frag
  $C/run_bscampp.sh $D fix "5000 10000" frag
  $C/run_bscampp.sh $D stock "5000 10000" frag
  $C/run_bscampp.sh $D fix 2000 frag
  $C/run_bscampp.sh $D "bug2only bug1only" 5000 frag
done
echo QUEUEDONE
