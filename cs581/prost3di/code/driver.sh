#!/bin/bash
# Run all methods on every BAliBASE set whose 3Di prediction is ready; loop until all sets are done.
CODE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-/opt/p3d/runs}
SETS=$(ls /opt/bb3/bb3_release/RV11/BB*.tfa /opt/bb3/bb3_release/RV12/BB*.tfa | xargs -n1 basename | sed 's/.tfa//')
while true; do
  left=0
  for s in $SETS; do
    if [ -s /opt/p3d/bb3/$s.3di.fa ]; then "$CODE/run_methods.sh" $s "$OUT"; else left=$((left+1)); fi
  done
  [ $left = 0 ] && break
  sleep 60
done
echo ALLDONE
