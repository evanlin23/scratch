#!/bin/bash
# De novo check on small datasets (100 sequences, 10 long), where MAFFT L-INS-i is
# affordable: does a de novo method recover the second domain / get misled by flanks?
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC=/opt/data/Datasets/ROSE/1000M4
for r in ${REPS:-0 1}; do
  for cond in "rand 0" "rand 1" "dom 0"; do
    set -- $cond
    ds=/opt/runs/lenhet/small_$1_m$2/R$r
    [ -f $ds/true.fasta ] || python3 $HERE/make_long.py $SRC/R$r/rose.aln.true.fasta $ds --n 100 --type $1 \
        --frac 0.1 --mult $2 --dom-src $SRC/R$(( (r + 10) % 20 ))/rose.aln.true.fasta --seed $((r + 1))
    for m in ${METHODS:-mafft upp upp-tfa emma mafft-add}; do echo "$ds $m"; done
  done
done | xargs -P ${P:-2} -L 1 bash -c 'python3 '"$HERE"'/run.py $0 $1 --threads 2 > /dev/null 2>> /opt/runs/lenhet/errors.log || echo "FAIL $0 $1" >> /opt/runs/lenhet/errors.log'
echo SMALL-DONE
