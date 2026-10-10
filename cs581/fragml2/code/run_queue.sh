#!/bin/bash
# Run queue.cmds 4 at a time (1 thread each). Restartable: finished arms are skipped by pipe.py.
cd "$(dirname "$0")"
export MLDATA=/opt/data/mlcache
xargs -P 4 -I CMD bash -c CMD < queue.cmds >> ../results/queue.log 2>&1
