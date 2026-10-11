#!/bin/bash
# Restartable driver: rerun after an interruption; finished subsets / merges are skipped.
# Proteins first (primary analysis). MUSCLE5 dropped on nucleotides: ~400-500 CPU s per 1,000-bp subset.
cd /home/user/scratch
export PYTHONPATH=cs581/code
B=cs581/basemeth; W=/opt/work/basemeth; PY="python3 $B/code/basemeth.py"
PROT=linsi,ginsi,muscle5,famsa,probcons,clustalo,kalign
DNA=linsi,ginsi,famsa,clustalo,kalign,prank
for round in 1 2 3 4 5 6 7 8 9 10; do
  $PY align $B/jobs.txt $W --methods $PROT
  $PY merge $B/jobs.txt $W $B/results/merge.jsonl --methods orig,$PROT
  $PY align $B/jobs_dna.txt $W --methods $DNA
  $PY merge $B/jobs_dna.txt $W $B/results/merge.jsonl --methods orig,$DNA
  if ! pgrep -f "basemeth.py prep" > /dev/null; then
    [ "$last" = 1 ] && break
    last=1
  fi
  sleep 60
done
echo ALLDONE
