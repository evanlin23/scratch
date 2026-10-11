#!/bin/bash
# wait_then.sh PID CMD...: run CMD after process PID exits
P=$1; shift; while kill -0 $P 2>/dev/null; do sleep 10; done; exec "$@"
