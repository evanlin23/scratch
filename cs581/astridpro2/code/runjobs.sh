#!/bin/bash
# Usage: runjobs.sh JOBFILE [PARALLEL] [METHODS]  -- restartable (bench.py skips finished (key, method)).
# Each line of JOBFILE is a list of bench.py arguments (shell-quoted).
H=$(cd "$(dirname "$0")" && pwd)
P=${2:-4}; M=${3:-}
export BENCH_TMP=/opt/tmp; mkdir -p $BENCH_TMP
EXTRA=""; [ -n "$M" ] && EXTRA="--methods $M"
export H EXTRA
tr '\n' '\0' < "$1" | xargs -0 -P "$P" -n 1 bash -c 'eval "/opt/mm/root/envs/gdl/bin/python $H/bench.py $1 $EXTRA" > /dev/null 2>&1' _
