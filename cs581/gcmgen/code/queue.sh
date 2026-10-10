#!/bin/bash
# Restartable: each line "REP VARIANT [VARIANT ...]" or "REP agree TOOL [...]"; finished rows are skipped.
# Usage: queue.sh LANEFILE   (reps under /opt/work/gcmgen/reps)
G=/home/user/scratch/cs581/gcmgen/code
cd /home/user/scratch/cs581/code
while read -r rep cmd rest; do
  [ -z "$rep" ] && continue; [[ $rep == \#* ]] && continue
  if [ "$cmd" = agree ]; then python3 $G/gg.py agree /opt/work/gcmgen/reps/$rep $rest
  else python3 $G/gg.py run /opt/work/gcmgen/reps/$rep $cmd $rest; fi
  bash $G/collect.sh
done < "$1"
echo QUEUE_DONE "$1"
