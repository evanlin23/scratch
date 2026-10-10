#!/bin/bash
# BSCAMPP's epa-ng: runs the EPA-ng binary named by $EPA_BIN (stock/fix/bug1only/bug2only builds).
# Logs one line per call (#tips via -t, wall time) to $EPA_CALLLOG if set.
s=$(date +%s.%N)
"${EPA_BIN:?set EPA_BIN}" "$@"; rc=$?
if [ -n "$EPA_CALLLOG" ]; then
  for ((i=1;i<=$#;i++)); do [ "${!i}" = -t ] && { j=$((i+1)); t=${!j}; }; done
  n=$(grep -o ',' "$t" | wc -l)
  echo -e "$((n+1))\t$(echo "$(date +%s.%N) - $s" | bc)\t$rc" >> "$EPA_CALLLOG"
fi
exit $rc
