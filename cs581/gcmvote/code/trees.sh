#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# FastTree -lg -gamma (protbench trees.py) on the true alignment and the vote.py outputs of a SIMHIGH replicate;
# nRF vs the true tree -> REP/trees.jsonl. Restartable (done rows skipped).
#   bash trees.sh NAME VARIANT_TAG [...]    e.g. bash trees.sh SIMHIGH_R1 magus es4 hard hard+mask:masked
# A tag "X:masked" uses vote/X/out.masked.fasta.
W=/opt/work/gcmvote; name=$1; shift
rep=$W/reps/$name; T=$W/trees/${name}_d0; mkdir -p $T
lvl=${name%_R*}; r=${name#*_R}
ln -sf $rep/true.fasta $T/true.fasta
ms="true"; [ -n "${NOTRUE:-}" ] && ms=""
for v in "$@"; do
  if [[ $v == *:masked ]]; then b=${v%:masked}; src=$rep/vote/$b/out.masked.fasta; m=${b//+/_}_masked
  else src=$rep/vote/$v/out.fasta; m=${v//+/_}; fi
  [ -s $src ] || { echo "missing $src"; continue; }
  ln -sf $src $T/$m.fasta; ms="$ms $m"
done
nice -n 5 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T \
  /opt/data/sim/$lvl/R$r/tree.nwk $rep/trees.jsonl $ms
cat $rep/trees.jsonl
