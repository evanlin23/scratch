#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh NAME : merge-only inputs of a replicate -> cs581/gcmvote/bank/NAME.tar.gz + MANIFEST.tsv row
set -eu
n=$1; lvl=${n%_R*}; r=${n#*_R}
R=/opt/work/gcmtrees/reps/$n; B=/home/user/scratch/cs581/gcmvote/bank; M=$B/MANIFEST.tsv
st=$(mktemp -d); mkdir $st/$n
cp -rL $R/inputs $R/sets $R/aligned $R/true.fasta $R/unaligned.fasta $R/magus.json $st/$n/
cp /opt/data/sim/$lvl/R$r/tree.nwk $st/$n/true_tree.nwk
[ -e $R/subsets.json ] && cp $R/subsets.json $st/$n/
mkdir -p $B; tar czf $B/$n.tar.gz -C $st $n; rm -rf $st
[ -f $M ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $M
grep -v "^$n	" $M > $M.tmp || true; mv $M.tmp $M
nseq=$(grep -c '>' $R/true.fasta)
printf '%s\tprotein\t%s\t10\t200\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..10}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..10}.fa)\n' \
  $n $nseq $(stat -c %s $B/$n.tar.gz) $(sha256sum $B/$n.tar.gz | cut -d' ' -f1) $n $n $n >> $M
