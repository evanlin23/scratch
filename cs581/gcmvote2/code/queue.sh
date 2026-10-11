#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Run run2.py on every rep in list file $1 (one path per line), 3 reps at a time; variants = remaining args.
# Restartable (run2.py skips variants already scored).
L=$1; shift
xargs -a "$L" -P 3 -I{} python3 /home/user/scratch/cs581/gcmvote2/code/run2.py {} "$@" > /dev/null
echo QUEUE_DONE "$L"
