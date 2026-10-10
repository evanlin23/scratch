#!/bin/bash
# Reproduce WITCH-NG (doi:10.1093/bioadv/vbad024) Table 2, 1000M3-HF: UPP and WITCH
# on all sequences, backbone = the 500 full-length sequences, queries = fragments.
# Two backbones: MAGUS-estimated (paper default) and true.
HERE="$(cd "$(dirname "$0")" && pwd)"
for r in ${REPS:-0 1 2}; do
  src=/opt/data/rosehf/high_frag/1000M3/R$r
  for bb in magus true; do
    ds=/opt/runs/lenhet/valid_1000M3HF_$bb/R$r
    mkdir -p $ds
    if [ ! -f $ds/backbone.txt ]; then
      cp $src/unaligned_all.txt $ds/unaligned.fasta
      cp $src/true_align_fragged.txt $ds/true.fasta
      grep ">" $src/unaligned_full.txt | sed 's/>//' > $ds/backbone.txt
      grep ">" $src/unaligned_frag.txt | sed 's/>//' > $ds/long.txt   # "long" = the fragments here
    fi
    for m in upp witch; do
      BACKBONE=$bb THREADS=${T:-4} python3 $HERE/run.py $ds $m --threads ${T:-4} > /dev/null 2>> /opt/runs/lenhet/valid_errors.log
    done
  done
done
echo VALID-DONE
