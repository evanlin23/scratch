#!/bin/bash
# After queue3: retry R0 stock b=9000 fragments (first attempt ran out of memory).
while pgrep -f "queue[23].sh" > /dev/null; do sleep 20; done
/home/user/scratch/cs581/epang/code/run_bscampp.sh /tmp/claude-0/ep/R0 "stock" "9000" "frag"
