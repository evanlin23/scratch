#!/bin/bash
# Whole-dataset baselines (restartable). Run after run_all.sh to limit CPU contention.
cd /home/user/scratch
export PYTHONPATH=cs581/code
B=cs581/basemeth; W=/opt/work/basemeth
while pgrep -f "basemeth.py (align|prep)" > /dev/null; do sleep 30; done
cat $B/jobs.txt $B/jobs_dna.txt > $W/jobs_all.txt
PY="python3 $B/code/baselines.py $W/jobs_all.txt $W $B/results/baselines.jsonl --cap 1800"
$PY --methods famsa
$PY --methods regressive
$PY --methods muscle5 --only BBA0039,BBA0067,BBA0101,BBA0154,BBA0190,SIMMOD_R1,SIMHIGH_R1,HF_rrm,HF_zf-CCHH,1000M2,16S.M
echo ALLDONE
