#!/bin/bash
# Dev round 2 (confidence-aware constraint), after dev round 1.
cd "$(dirname "$0")/.."
while pgrep -x xargs > /dev/null; do sleep 60; done
xargs -P 4 -I CMD bash -c CMD < code/dev2.cmds
