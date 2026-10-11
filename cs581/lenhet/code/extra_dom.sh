#!/bin/bash
# Extra replicates of the second-domain condition (more power for the paired test).
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC=/opt/data/Datasets/ROSE/1000M4
for r in ${REPS:-3 4 5}; do
  ds=/opt/runs/lenhet/dom_f0.1_m0/R$r
  [ -f $ds/true.fasta ] || python3 $HERE/make_long.py $SRC/R$r/rose.aln.true.fasta $ds --n 500 --type dom \
      --frac 0.1 --dom-src $SRC/R$(( (r + 10) % 20 ))/rose.aln.true.fasta --seed $((r + 1))
  for m in ${METHODS:-upp emma mafft-add upp-tfa emma-tfa witch}; do echo "$ds $m"; done
done | xargs -P ${P:-2} -L 1 bash -c 'python3 '"$HERE"'/run.py $0 $1 --threads 1 > /dev/null 2>> /opt/runs/lenhet/errors.log || echo "FAIL $0 $1" >> /opt/runs/lenhet/errors.log'
echo EXTRA-DONE
