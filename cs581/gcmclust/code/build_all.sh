#!/bin/bash
# Lane: build graph + control merge for each rep in the list file (restartable).
while read -r n; do python3 /home/user/scratch/cs581/gcmclust/code/gc.py build /opt/work/gcmclust/reps/$n || echo BUILD_FAIL $n; done < "$1"
echo BUILD_DONE "$1"
