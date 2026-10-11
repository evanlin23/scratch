#!/bin/bash
# Usage: sim_batch.sh SHAPE INTERNAL_MEAN STMODE REP_FROM REP_TO PARALLEL
cd "$(dirname "$0")"
S=$1; IM=$2; ST=$3; A=$4; B=$5; P=${6:-3}
OUT=/opt/gtmdata/simq
seq $A $B | xargs -P $P -I{} bash -c "QUICK=1 INTERNAL_MEAN=$IM python3 sim.py $S {} $OUT 200 50 $ST > /dev/null 2>&1; python3 sim_ml.py $OUT/${S}_n200_m50_i$IM/{} $ST 4 > /dev/null 2>&1; echo done $S $IM $ST {}"
