#!/bin/bash
# Run bbe.py jobs from a lane file serially: each line "REP VARIANT [VARIANT ...]" (run, then diag).
# Usage: queue.sh LANEFILE
cd "$(dirname "$0")"
while read -r rep vars; do
  [ -z "$rep" ] && continue
  python3 bbe.py run /home/user/work/reps/$rep $vars >> /home/user/work/queue_$(basename $1).log 2>&1
  python3 bbe.py diag /home/user/work/reps/$rep $vars >> /home/user/work/queue_$(basename $1).log 2>&1
done < "$1"
