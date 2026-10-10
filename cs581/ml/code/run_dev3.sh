#!/bin/bash
# Dev round 3 (cheap polish, stricter backbone threshold), after round 2.
cd "$(dirname "$0")/.."
while pgrep -f "code/run_dev2.sh" > /dev/null; do sleep 60; done
xargs -P 4 -I CMD bash -c CMD < code/dev3.cmds
