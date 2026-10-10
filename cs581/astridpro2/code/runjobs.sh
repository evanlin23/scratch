#!/bin/bash
# Usage: runjobs.sh JOBFILE [PARALLEL] [METHODS]  -- restartable (bench.py skips finished (key, method)).
H=$(cd "$(dirname "$0")" && pwd)
P=${2:-4}; M=${3:-}
export BENCH_TMP=/opt/tmp; mkdir -p $BENCH_TMP
EXTRA=""; [ -n "$M" ] && EXTRA="--methods $M"
xargs -P "$P" -I{} bash -c "/opt/mm/root/envs/gdl/bin/python $H/bench.py {} $EXTRA > /dev/null 2>&1" < "$1"
