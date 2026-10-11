#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh NAME : merge-only inputs of a gcmtrees rep -> cs581/gcmvote/bank/NAME.tar.gz + MANIFEST.tsv row
set -eu
n=$1; G=/opt/work/gcmtrees/reps/$n; r=${n#*_R}; B=/home/user/scratch/cs581/gcmvote/bank
S=$(mktemp -d); mkdir $S/$n
cp -r $G/inputs $G/sets $G/aligned $G/true.fasta $G/unaligned.fasta $G/magus.json $S/$n/
cp /opt/data/sim/SIMHIGH/R$r/tree.nwk $S/$n/true_tree.nwk
cp /opt/work/gcmvote/reps/$n/vote/magus/subsets.json $S/$n/subsets.json
mkdir -p $B; tar -C $S -czf $B/$n.tar.gz $n; rm -rf $S
M=$B/MANIFEST.tsv
[ -f $M ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $M
sed -i "/^$n\t/d" $M
nb=$(ls $G/inputs/backbones | wc -l); ns=$(grep -c '>' $G/unaligned.fasta)
bs=$(grep -c '>' $G/inputs/backbones/backbone_1_mafft.txt)
printf '%s\tprotein\t%s\t%s\t%s\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..%s}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..%s}.fa)\n' \
  $n $ns $nb $bs $(stat -c %s $B/$n.tar.gz) $(sha256sum $B/$n.tar.gz | cut -d' ' -f1) $n $n $nb $n $nb >> $M
tail -1 $M
