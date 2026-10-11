#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Rep bank: merge-only inputs of each h4c replicate -> bank/NAME.tar.gz (+ MANIFEST.tsv line). Merged outputs excluded.
W=/opt/work/gcmvote; B=/home/user/scratch/cs581/gcmvote/bank; mkdir -p $B
[ -f $B/MANIFEST.tsv ] || printf "name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\tref_path\tbackbones_path\n" > $B/MANIFEST.tsv
for name in "$@"; do
  rep=$W/reps/$name; t=$(mktemp -d); d=$t/$name; mkdir -p $d
  cp -r $rep/inputs $rep/sets $rep/aligned $rep/true.fasta $rep/unaligned.fasta $rep/magus.json $d/
  src=/opt/data/Datasets/ROSE/$name/R0
  cp $src/rose.tt $d/true.tree; cp $src/random.tree $d/random.tree
  tar -C $t -czf $B/$name.tar.gz $name
  nseq=$(grep -c '>' $d/unaligned.fasta); nbb=$(ls $d/inputs/backbones | wc -l)
  bsz=$(grep -c '>' $d/inputs/backbones/backbone_1_mafft.txt)
  bytes=$(stat -c %s $B/$name.tar.gz); sha=$(sha256sum $B/$name.tar.gz | cut -d' ' -f1)
  printf "%s\tdna\t%s\t%s\t%s\tyes\t%s\t%s\t%s\t%s\n" $name $nseq $nbb $bsz $bytes $sha $name/true.fasta $name/inputs/backbones/ >> $B/MANIFEST.tsv
  rm -rf $t
done
