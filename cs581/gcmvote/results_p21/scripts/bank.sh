#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh SIMHIGH_Rk : merge-only inputs tarball + MANIFEST row
set -eu
n=$1; r=${n#*_R}; lvl=${n%_R*}
R=/opt/work/gcmtrees/reps/$n; B=/home/user/scratch/cs581/gcmvote/bank; S=/opt/work/bankstage
rm -rf $S/$n; mkdir -p $S/$n $B
cp -rL $R/inputs $R/sets $R/aligned $R/true.fasta $R/unaligned.fasta $R/magus.json $S/$n/
cp /opt/data/sim/$lvl/R$r/tree.nwk $S/$n/true_tree.nwk
cp /opt/work/gcmvote/reps/$n/vote/magus/subsets.json $S/$n/subsets.json
tar czf $B/$n.tar.gz -C $S $n
[ -f $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
nseq=$(grep -c '>' $S/$n/true.fasta); bs=$(grep -c '>' $S/$n/inputs/backbones/backbone_1_mafft.txt)
printf '%s\tprotein\t%s\t10\t%s\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..10}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..10}.fa)\n' \
  $n $nseq $bs $(stat -c %s $B/$n.tar.gz) $(sha256sum $B/$n.tar.gz | cut -c1-64) $n $n $n >> $B/MANIFEST.tsv
