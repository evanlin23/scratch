#!/bin/bash
# Extra replicates for the end-to-end BSCAMPP comparison (after queues 2-4 and replicate prep).
C=/home/user/scratch/cs581/epang/code; E=/tmp/claude-0/ep
while pgrep -f "queue[234].sh" > /dev/null; do sleep 20; done
for r in R2 R3; do
  until grep -q REPDONE $E/$r/rx.out 2>/dev/null; do sleep 20; done
  $C/run_bscampp.sh $E/$r "stock fix" "2000 5000" "frag full"
  $C/run_bscampp.sh $E/$r "fix" "9000" "frag"
done
