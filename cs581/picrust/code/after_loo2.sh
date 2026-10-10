#!/bin/bash
# wait for after_loo.sh, then run the V3-130nt leave-out query set (s2) with stock and fix, and score it
while pgrep -f after_loo.sh >/dev/null; do sleep 15; done
C=/home/user/scratch/cs581/picrust/code
for v in stock fix; do $C/run_loo.sh /opt/work/loo/s2 $v; done
/opt/mm/root/envs/picrust2/bin/python $C/score_loo.py /opt/work/loo/s2 > /opt/work/loo/s2/summary.json 2> /opt/work/loo/s2/score.err
