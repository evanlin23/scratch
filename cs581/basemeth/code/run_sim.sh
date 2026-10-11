#!/bin/bash
# SIM protein sets first (orchestrator priority), then the general driver.
cd /home/user/scratch
export PYTHONPATH=cs581/code
B=cs581/basemeth; W=/opt/work/basemeth; PY="python3 $B/code/basemeth.py"
PROT=linsi,ginsi,muscle5,famsa,probcons,clustalo,kalign
for round in 1 2 3 4 5 6; do
  $PY align $B/jobs_sim.txt $W --methods $PROT
  $PY merge $B/jobs_sim.txt $W $B/results/merge.jsonl --methods orig,$PROT
  n=$(ls -d $W/SIM*_R*/prep.json 2>/dev/null | wc -l); [ $n -ge 4 ] && [ $round -ge 2 ] && break
  sleep 60
done
cs581/basemeth/code/run_all.sh
