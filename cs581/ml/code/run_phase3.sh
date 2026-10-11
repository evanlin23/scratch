#!/bin/bash
# Phase 3 (expensive RAxML-NG pilots): one single-threaded job per core, 4 at a time.
cd "$(dirname "$0")/.."
while pgrep -f "code/run_phase2.sh" > /dev/null || pgrep -f "code/run_probe.sh" > /dev/null; do sleep 30; done
xargs -P 4 -I CMD bash -c CMD < code/phase3.cmds
