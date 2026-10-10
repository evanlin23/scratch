#!/bin/bash
# Frozen methods (FastTree backbone, support cutoff 0.8, constrained RAxML-NG, + fast polish) on the
# held-out test sets, plus the estimated-alignment (UPP) arm; after the test baselines.
cd "$(dirname "$0")/.."
while pgrep -f "code/run_test_base.sh" > /dev/null; do sleep 60; done
xargs -P 4 -I CMD bash -c CMD < code/test_final.cmds
