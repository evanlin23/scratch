#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# start each queued tree when fewer than 4 FastTree processes run; skip methods already in trees.jsonl
J=/opt/work/gcmtrees/reps/SIMHIGH_R20/trees.jsonl
while read m src; do
  grep -q "\"method\": \"$m\"" $J 2>/dev/null && continue
  pgrep -f "trees.py .* $m\$" >/dev/null && continue
  while [ $(pgrep -c -x FastTree) -ge 4 ]; do sleep 15; done
  nohup bash /opt/work/tree1.sh $m $src >/dev/null 2>&1 &
  sleep 5
done < /opt/work/queue.txt
wait; echo SCHED_DONE
