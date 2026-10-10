#!/bin/bash
# Controls requested after the diagnosis (R0 nested): bug-1-only / bug-2-only builds, fixed build with
# forced rate scalers and with --no-heur, and complete --rate-scalers off at >2000 tips.
C=/home/user/scratch/cs581/epang/code; cd /tmp/claude-0/ep/R0
while pgrep -f "queue[2345].sh" > /dev/null; do sleep 20; done
export EPA_FIX=/opt/tools/epa/epa-ng-fix2
KS="500 1000 2000" $C/run_nested.sh . "fix_rson" "frag" 2
KS="3000 5000" $C/run_nested.sh . "shift simd rs_off" "frag" 2
KS="2000 3000 5000" $C/run_nested.sh . "fix_noheur" "frag" 2
