#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Rep bank: merge-only inputs of each h6b rep dir -> bank/<REP>.tar.gz (+ MANIFEST.tsv line).
# Contents: inputs/ (subalignments + 20 L-INS-i backbones), aligned/ (gg.py's backbone copies), bb5/bb10/bb20,
# sets/ (unaligned backbone sets), true.fasta (reference), unaligned.fasta, magus.json, subsets.json
# (MAGUS's subset order from vote/hard@B5), true_tree.nwk (AliSim). No merged outputs.
set -eu
W=/opt/work/h6b/reps
BANK=$(cd "$(dirname "$0")/../bank" && pwd)
M=$BANK/MANIFEST.tsv
[ -f $M ] || printf 'name\ttype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $M
for rep in "$@"; do
  d=$W/$rep; lvl=${rep%_R*}; r=${rep#*_R}
  st=$(mktemp -d); mkdir $st/$rep
  cp -r $d/inputs $d/aligned $d/bb5 $d/bb10 $d/bb20 $d/sets $d/true.fasta $d/unaligned.fasta $d/magus.json $st/$rep/
  cp $d/vote/hard@B5/subsets.json $st/$rep/subsets.json
  cp /opt/data/sim/$lvl/R$r/tree.nwk $st/$rep/true_tree.nwk
  tar -C $st -czf $BANK/$rep.tar.gz $rep; rm -rf $st
  bytes=$(stat -c %s $BANK/$rep.tar.gz)
  if [ $bytes -gt 52428800 ]; then echo "SKIP $rep ($bytes bytes > 50 MB)"; rm $BANK/$rep.tar.gz; continue; fi
  nseq=$(grep -c '>' $d/unaligned.fasta); nbb=$(ls $d/inputs/backbones | wc -l)
  bsz=$(grep -c '>' $d/inputs/backbones/backbone_1_mafft.txt)
  printf '%s\tprotein\t%s\t%s\t%s\tyes\t%s\t%s\t%s\t%s\n' $rep $nseq $nbb $bsz $bytes $(sha256sum $BANK/$rep.tar.gz | cut -d' ' -f1) \
    $rep/true.fasta "$rep/inputs/backbones/backbone_{1..$nbb}_mafft.txt (B-subsets: $rep/bb{5,10,20}/)" >> $M
  echo "$rep $bytes"
done
