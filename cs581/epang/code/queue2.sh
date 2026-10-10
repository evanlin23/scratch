#!/bin/bash
# End-to-end BSCAMPP queue (run alone for clean runtimes), then remaining nested runs.
C=/home/user/scratch/cs581/epang/code
E=/tmp/claude-0/ep
$C/run_bscampp.sh $E/R0 "stock fix" "2000 5000" "frag full"
$C/run_bscampp.sh $E/R0 "stock fix" "9000 1000" "frag full"
cd $E/R0 && ONLY_ALL=1 THREADS=4 EPA_FIX=/opt/tools/epa/epa-ng-fix2 $C/run_nested.sh . "auto fix" "frag full" 1
$C/run_bscampp.sh $E/R1 "stock fix" "2000 5000 9000" "frag full"
cd $E/R0 && $C/run_nested.sh . "noheur" "frag" 2
