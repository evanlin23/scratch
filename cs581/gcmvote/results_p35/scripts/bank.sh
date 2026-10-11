#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh NAME : merge-only inputs tarball + MANIFEST row
n=$1; lvl=${n%_R*}; r=${n#*_R}; G=/opt/work/gcmtrees/reps/$n
B=/home/user/scratch/cs581/gcmvote/bank; mkdir -p $B
st=$(mktemp -d); mkdir $st/$n
cp -r $G/inputs $G/sets $G/aligned $G/true.fasta $G/unaligned.fasta $G/magus.json $st/$n/
[ -f $G/subsets.json ] && cp $G/subsets.json $st/$n/
cp /opt/data/sim/$lvl/R$r/tree.nwk $st/$n/true_tree.nwk
tar czf $B/$n.tar.gz -C $st $n; rm -rf $st
M=$B/MANIFEST.tsv
[ -f $M ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $M
nseq=$(grep -c '>' $G/unaligned.fasta); nbb=$(ls $G/inputs/backbones | wc -l)
bsz=$(grep -c '>' $G/inputs/backbones/backbone_1_mafft.txt)
printf '%s\tprotein\t%s\t%s\t%s\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..%s}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..%s}.fa)\n' \
  $n $nseq $nbb $bsz $(stat -c%s $B/$n.tar.gz) $(sha256sum $B/$n.tar.gz | cut -c1-64) $n $n $nbb $n $nbb >> $M
tail -1 $M
