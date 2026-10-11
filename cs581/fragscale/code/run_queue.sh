#!/bin/bash
# Run the lines of a queue file (pipe.py arguments) with P parallel workers. Restartable.
#   bash run_queue.sh QUEUE.cmds P [extra pipe.py args]
HERE="$(cd "$(dirname "$0")" && pwd)"
Q=$1; P=$2; shift 2
export MLDATA=${MLDATA:-/opt/data/fscache}
PY=/opt/mm/root/envs/fml/bin/python
grep -v '^#' "$Q" | grep . | xargs -P "$P" -I{} sh -c "$PY $HERE/pipe.py {} $*"
