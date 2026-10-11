#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Queue cut to 1000L3 + 1000M3: stop run_base.sh (and its children) once 1000M3's baselines are in.
until grep -q "BASE_DONE 1000M3_R0$" /opt/work/h4/base.log; do sleep 15; done
for p in $(pgrep -f run_base.sh); do pkill -TERM -P $p; kill $p; done
echo ALL_BASE_DONE >> /opt/work/h4/base.log
