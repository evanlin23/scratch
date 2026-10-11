#!/bin/bash
# Placed on PATH as "epa-ng" so PICRUSt2's place_seqs.py calls the chosen build.
# EPA_BIN selects the binary; each call's wall time and args are appended to $EPA_LOG.
s=$(date +%s.%N)
"$EPA_BIN" "$@"; rc=$?
e=$(date +%s.%N)
[ -n "$EPA_LOG" ] && echo -e "$(basename $(dirname $EPA_BIN))\t$(echo "$e - $s" | bc)\t$*" >> "$EPA_LOG"
exit $rc
