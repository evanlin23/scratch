#!/bin/bash
while pgrep -f "queue[23456].sh" > /dev/null; do sleep 20; done
/home/user/scratch/cs581/epang/code/timing_full.sh /tmp/claude-0/ep/R0
