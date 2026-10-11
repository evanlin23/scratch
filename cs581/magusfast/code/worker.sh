#!/bin/bash
# Restartable queue: runs each "STEP NAME" line of queue.txt in order (finished steps are skipped by
# mf_bench.py). Rerun after the session stops background commands.
cd "$(dirname "$0")"
while true; do
  line=$(grep -v '^#' queue.txt | grep -v -x -F -f done.txt 2>/dev/null | head -1)
  [ -z "$line" ] && { echo "queue empty"; break; }
  echo "$(date +%T) START $line"
  if python3 mf_bench.py $line >> worker.log 2>&1; then
    echo "$line" >> done.txt; echo "$(date +%T) DONE $line"
  else
    echo "$(date +%T) FAIL $line"; echo "$line" >> done.txt; echo "FAILED $line" >> failed.txt
  fi
done
