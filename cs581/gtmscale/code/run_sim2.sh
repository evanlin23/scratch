#!/bin/bash
# second pass after run_sim.sh: missed jobs + one iteration round (guide = previous Blend-FT tree)
cd "$(dirname "$0")"
export THREADS=${THREADS:-2}
while pgrep -x run_sim.sh >/dev/null || ps -ef | grep -v grep | grep -q "run_sim.sh$"; do sleep 30; done
S=/opt/gtms/sim
python3 pipe.py $S/n2000_i0.01/r1 kmer 500 >> $S/log2.txt 2>&1
for d in $S/n2000_i0.01/r* $S/n2000_i0.02/r* $S/n5000_i0.01/r*; do
  python3 pipe.py $d ftfast+it 500 gtm,blendft,polishft >> $S/log2.txt 2>&1
done
for r in 1 2 3 4; do python3 pipe.py $S/n2000_i0.01/r$r kmer+it 500 gtm,blendft,polishft >> $S/log2.txt 2>&1; done
echo DONE2 >> $S/log2.txt
