#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Rep bank: tar the merge-only inputs of a rep dir (no merged outputs) and append a MANIFEST.tsv line.
#   bash bank.sh NAME REPDIR DTYPE [TRUE_TREE]
name=$1; rep=$2; dtype=$3; tree=$4
B=/home/user/scratch/cs581/gcmvote/bank; tmp=$(mktemp -d); d=$tmp/$name; mkdir -p $d
cp -r $rep/inputs $rep/true.fasta $rep/unaligned.fasta $rep/magus.json $d/
[ -d $rep/sets ] && cp -r $rep/sets $d/
for x in $rep/bb5 $rep/bb20; do [ -d $x ] && cp -r $x $d/; done
has_tree=0; [ -n "$tree" ] && { cp $tree $d/true_tree.tre; has_tree=1; }
tar -C $tmp -czf $B/$name.tar.gz $name
bytes=$(stat -c %s $B/$name.tar.gz)
if [ $bytes -gt 52428800 ]; then echo "SKIP $name: $bytes bytes > 50 MB"; rm $B/$name.tar.gz; rm -rf $tmp; exit 0; fi
nseq=$(grep -c '>' $d/unaligned.fasta); nbb=$(ls $d/inputs/backbones | wc -l)
bbsize=$(grep -c '>' $d/inputs/backbones/backbone_1_mafft.txt)
sha=$(sha256sum $B/$name.tar.gz | cut -d' ' -f1)
[ -f $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $name $dtype $nseq $nbb $bbsize $has_tree $bytes $sha \
  $name/true.fasta $name/inputs/backbones/ >> $B/MANIFEST.tsv
rm -rf $tmp; echo "$name $bytes"
