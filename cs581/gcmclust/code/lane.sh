#!/bin/bash
# Lane: for each rep in LISTFILE (waits for its graph), run the variants VARIANTS... (restartable; per-variant locks
# let several lanes share replicates).
L=$1; shift
for n in $(cat $L); do
  until [ -f /opt/work/gcmclust/reps/$n/graph.npz ]; do sleep 30; done
  python3 /home/user/scratch/cs581/gcmclust/code/gc.py run /opt/work/gcmclust/reps/$n "$@"
done
echo LANE_DONE $L
