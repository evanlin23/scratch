#!/bin/bash
# Phase 4: the remaining 1000M1-HF replicates (n=5 as in Park et al. 2021), after phase 3.
cd "$(dirname "$0")/.."
while pgrep -x xargs > /dev/null; do sleep 60; done
xargs -P 4 -I CMD bash -c CMD < code/phase4.cmds
