#!/bin/bash
# simulation + RNASim10K jobs, sequential; restartable (pipe.py skips finished steps)
cd "$(dirname "$0")"
export THREADS=${THREADS:-2} IQ_CAP=${IQ_CAP:-3600}
S=/opt/gtms/sim
job() { # n imean rep guide maxsub
  d=$S/n$1_i$2/r$3
  python3 simgen.py $d $1 $2 $3 && python3 pipe.py $d $4 $5 >> $S/log_n$1_i$2_r$3_$4.txt 2>&1
}
for r in 1 2 3 4 5 6; do job 2000 0.01 $r ftfast 500; done
for r in 1 2 3 4; do job 2000 0.01 $r kmer 500; done
for r in 1 2 3 4; do job 2000 0.02 $r ftfast 500; done
for r in 1 2; do job 5000 0.01 $r ftfast 500; done
for r in 0 1; do
  d=/opt/gtms/rnasim10k/R$r; mkdir -p $d
  [ -f $d/aln.fa ] || ln -sf /opt/data/Datasets/RNASim/10000/R$r/true_align.txt $d/aln.fa
  [ -f $d/true.tre ] || ln -sf /opt/data/Datasets/RNASim/10000/R$r/true_tree.tre $d/true.tre
  python3 pipe.py $d ftfast 1000 full_ft,gtm,blendft,polishft >> $d/log.txt 2>&1
done
echo ALLDONE >> $S/done.txt
