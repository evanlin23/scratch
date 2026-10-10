#!/bin/bash
# PASTA 1.8.3 with its built-in subset aligners (restartable; finished rows are skipped).
cd /home/user/scratch
export PYTHONPATH=cs581/code
B=cs581/basemeth; W=/opt/work/basemeth
while pgrep -f "baselines.py.*pasta" > /dev/null; do sleep 30; done  # an earlier copy still running
cat $B/jobs.txt $B/jobs_dna.txt > $W/jobs_all.txt
PY="python3 $B/code/baselines.py $W/jobs_all.txt $W $B/results/pasta.jsonl --cap 3600"
$PY --iter 1 --methods pasta-mafft,pasta-probcons,pasta-prank --only BBA0039,1000M2
$PY --iter 3 --methods pasta-mafft --only BBA0039,1000M2
$PY --iter 1 --methods pasta-mafft,pasta-probcons,pasta-prank --only SIMMOD_R1
$PY --iter 3 --methods pasta-mafft --only SIMMOD_R1
echo ALLDONE
