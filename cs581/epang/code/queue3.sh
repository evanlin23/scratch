#!/bin/bash
# Waits for queue2, then re-runs the R0 b=5000 BSCAMPP runs (first attempt ran out of memory).
while pgrep -f queue2.sh > /dev/null; do sleep 20; done
/home/user/scratch/cs581/epang/code/run_bscampp.sh /tmp/claude-0/ep/R0 "stock fix" "5000" "frag full"
