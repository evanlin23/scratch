#!/bin/bash
# Test-set baselines (independent of the method frozen on dev), after dev round 3.
cd "$(dirname "$0")/.."
while pgrep -f "code/run_dev3.sh" > /dev/null; do sleep 60; done
xargs -P 4 -I CMD bash -c CMD < code/test_base.cmds
