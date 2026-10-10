#!/bin/bash
# RNASim 50K (Recursive-MAGUS release, IDB-1048258, RNASim/50000_R_1/R0): 1000 queries,
# FastTree-2 backbone (GTR+gamma). A RAxML-NG --evaluate on the whole 49K backbone does not fit in
# 15 GB, so GTR+G parameters come from RAxML-NG --evaluate on a 5,000-leaf pruned subtree
# (model_from_subtree.py) and branch lengths stay FastTree's.
# Usage: prep_rnasim.sh <true_align.txt> <true_tree.tre> <outdir> [threads=1]
A=$1; TT=$2; O=$3; T=${4:-1}
P=/opt/mm/root/envs/place/bin; C=$(cd "$(dirname "$0")" && pwd)
$P/python $C/prep.py $A $O 1000 1 > $O.prep.log
cd $O
tr -d ' \n\r\t' < $TT > true_clean.tre; echo >> true_clean.tre; echo $O/true_clean.tre > true_tree.txt
OMP_NUM_THREADS=$T /usr/bin/time -v -o fasttree.time $P/FastTreeMP -nt -gtr -gamma backbone.fa > fasttree.tre 2> fasttree.log
$P/python $C/model_from_subtree.py . 5000 $T > rx.out 2>&1
echo PREPDONE >> rx.out
