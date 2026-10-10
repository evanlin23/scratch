#!/bin/bash
# Build the benchmark grid and run every (dataset, method) job; restartable
# (run.py skips jobs that already have score.json).
#   bash grid.sh [JOBS_IN_PARALLEL] [THREADS_PER_JOB]
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
P=${1:-2}; T=${2:-2}
ROOT=/opt/runs/lenhet
SRC=/opt/data/Datasets/ROSE/1000M4   # 1000M2 is near saturation (median p-distance 0.69)
REPS="${REPS:-0 1 2 3}"
METHODS="${METHODS:-upp witch emma mafft-add mafft-addlong upp-tfa emma-tfa}"
mkdir -p $ROOT
jobs=$ROOT/jobs.txt; : > $jobs
for r in $REPS; do
  aln=$SRC/R$r/rose.aln.true.fasta
  other=$SRC/R$(( (r + 10) % 20 ))/rose.aln.true.fasta   # family B for "dom"
  for cond in "rand f0.1 m0" "rand f0.1 m0.5" "rand f0.1 m1" "rand f0.1 m3" "dom f0.1 m0"; do
    set -- $cond
    ds=$ROOT/$1_$2_$3/R$r
    if [ ! -f $ds/true.fasta ]; then
      python3 $HERE/make_long.py $aln $ds --n 500 --type $1 --frac ${2#f} --mult ${3#m} \
        --dom-src $other --seed $((r + 1))
    fi
    for m in $METHODS; do echo "$ds $m" >> $jobs; done
  done
done
# de novo methods on the main conditions only (MAGUS is the slow one)
for r in ${DREPS:-0 1 2}; do for c in rand_f0.1_m0 rand_f0.1_m1 dom_f0.1_m0; do
  for m in ${DENOVO:-mafft magus}; do echo "$ROOT/$c/R$r $m" >> $jobs; done; done; done
cat $jobs | xargs -P $P -L 1 bash -c 'python3 '"$HERE"'/run.py $0 $1 --threads '"$T"' > /dev/null 2>> '"$ROOT"'/errors.log || echo "FAIL $0 $1" >> '"$ROOT"'/errors.log; echo "done $0 $1" >> '"$ROOT"'/progress.log'
echo GRID-DONE >> $ROOT/progress.log
