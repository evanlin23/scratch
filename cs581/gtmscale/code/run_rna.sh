#!/bin/bash
# RNASim 10K (true alignment, true tree), restartable
cd "$(dirname "$0")"
export THREADS=${THREADS:-2} IQ_CAP=${IQ_CAP:-5400}
for r in ${REPS:-0 1}; do
  d=/opt/gtms/rnasim10k/R$r; mkdir -p $d
  [ -e $d/aln.fa ] || ln -sf /opt/data/Datasets/RNASim/10000/R$r/true_align.txt $d/aln.fa
  [ -e $d/true.tre ] || ln -sf /opt/data/Datasets/RNASim/10000/R$r/true_tree.tre $d/true.tre
  python3 pipe.py $d ftfast 1000 full_ft,gtm,blendft,polishft >> $d/log.txt 2>&1
  python3 pipe.py $d ftfast+it 1000 gtm,blendft >> $d/log.txt 2>&1
done
# full IQ-TREE run separately (iq_full.py), alone, because it needs >13 GB
