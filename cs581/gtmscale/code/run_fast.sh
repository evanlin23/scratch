#!/bin/bash
cd "$(dirname "$0")"
for d in /opt/gtms/sim/n*_i*/r*; do
  for g in ftfast kmer ftfast+it; do [ -d $d/${g}_m500 ] && python3 pipe.py $d $g 500 blendfast >> /opt/gtms/sim/log_fast.txt 2>&1; done
done
echo DONEFAST >> /opt/gtms/sim/log_fast.txt
