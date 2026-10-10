#!/bin/bash
# Why-not lane, one FastTree at a time per dataset: zero-SPFP refinements (split_*) and FastTree search noise
# (same alignment, sequences shuffled with seed 1: *_shuf).   bash extra.sh DATASET [...]
W=/opt/work/gcmtrees; G=/home/user/scratch/cs581/gcmtrees
for n in "$@"; do
  T=$W/trees/${n}_d0; rep=$W/reps/$n; lvl=${n%_R*}; r=${n#*_R}
  until [ -s $T/recipe.fasta.fasttree.nwk ] || [ -s $T/recipe.fasttree.nwk ]; do sleep 60; done
  [ -s $T/split_magus.fasta ] || python3 $G/code/refine.py $T/true.fasta $T/magus.fasta $T/split_magus.fasta
  [ -s $T/split_recipe.fasta ] || python3 $G/code/refine.py $T/true.fasta $T/recipe.fasta $T/split_recipe.fasta
  for m in magus recipe; do [ -s $T/${m}_shuf.fasta ] || python3 $G/code/shuffle.py $T/$m.fasta $T/${m}_shuf.fasta 1; done
  nice -n 10 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T \
    /opt/data/sim/$lvl/R$r/tree.nwk $rep/trees.jsonl ${EXTRA_M:-split_magus split_recipe magus_shuf recipe_shuf} >/dev/null
  bash $G/code/collect.sh
done
echo EXTRA_DONE
