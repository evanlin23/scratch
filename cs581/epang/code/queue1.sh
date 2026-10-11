#!/bin/bash
# Nested-experiment queue used for the pilot (R0).
D=$1
export EPA_FIX=/opt/tools/epa/epa-ng-fix2
C=/home/user/scratch/cs581/epang/code
$C/run_nested.sh $D "fix auto" "full" 2
$C/run_nested.sh $D "nopremask noheur" "frag" 2
