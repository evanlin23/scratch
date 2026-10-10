#!/bin/bash
# Prepare one RNASim 10K replicate: split, FastTree backbone, RAxML-NG model/branch lengths.
# Usage: prep_rep.sh <rep, e.g. R1> <workdir> [threads=1]
R=$1; W=$2; T=${3:-1}
P=/opt/mm/root/envs/place/bin
cd $W
$P/python /home/user/scratch/cs581/epang/code/prep.py /opt/data/Datasets/RNASim/10000/$R/true_align.txt $R > $R.prep.log
cd $R
echo /opt/data/Datasets/RNASim/10000/$R/true_tree.tre > true_tree.txt
OMP_NUM_THREADS=$T $P/FastTreeMP -nt -gtr -gamma backbone.fa > fasttree.tre 2> fasttree.log
$P/raxml-ng --evaluate --msa backbone.fa --tree fasttree.tre --model GTR+G --prefix rx --threads $T --redo > rx.out 2>&1
echo REPDONE >> rx.out
