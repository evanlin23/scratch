#!/bin/bash
# copy small per-instance files (flagged sets + rogue lists + instance info) into results/detect/
S=/tmp/claude-0/-home-user-scratch/690a2cdc-90a2-5a40-b34f-c88f48d84cea/scratchpad
D="$(dirname "$0")/../results/detect"; mkdir -p "$D"
for w in $S/work/*/; do
  i=$(basename $w)
  [ -f $w/detect.json ] || continue
  python3 -c "
import json; j=json.load(open('$w/detect.json')); j['rogues']=open('$S/inst/$i/rogues.txt').read().split()
j['info']=json.load(open('$S/inst/$i/info.json')); j['info'].pop('anchors',None)
json.dump(j,open('$D/$i.json','w'))"
done
