#!/bin/bash
# Why-not diagnostic: FastTree on the zero-SPFP refinements (refine.py) of MAGUS's and the recipe's alignments.
#   bash diag.sh DATASET [...]     (after trees.sh has made trees/<DATASET>_d0)
W=/opt/work/gcmtrees; G=/home/user/scratch/cs581/gcmtrees
for n in "$@"; do
  T=$W/trees/${n}_d0; rep=$W/reps/$n; lvl=${n%_R*}; r=${n#*_R}
  [ -s $T/split_magus.fasta ] || python3 $G/code/refine.py $T/true.fasta $T/magus.fasta $T/split_magus.fasta
  [ -s $T/split_recipe.fasta ] || python3 $G/code/refine.py $T/true.fasta $T/recipe.fasta $T/split_recipe.fasta
  nice -n 10 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T \
    /opt/data/sim/$lvl/R$r/tree.nwk $rep/trees.jsonl split_magus split_recipe >/dev/null
  bash $G/code/collect.sh
done
echo DIAG_DONE
