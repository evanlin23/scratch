#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh SIMHIGH_Rk : merge-only inputs tarball + MANIFEST row
set -eu
n=$1; lvl=${n%_R*}; r=${n#*_R}
R=/opt/work/gcmtrees/reps/$n; B=/home/user/scratch/cs581/gcmvote/bank
ST=/opt/work/bankstage; rm -rf $ST/$n; mkdir -p $ST/$n $B
cp -r $R/inputs $R/sets $R/aligned $R/true.fasta $R/unaligned.fasta $R/magus.json $ST/$n/
[ -f $R/subsets.json ] && cp $R/subsets.json $ST/$n/
cp /opt/data/sim/$lvl/R$r/tree.nwk $ST/$n/true_tree.nwk
tar czf $B/$n.tar.gz -C $ST $n
[ -f $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
nseq=$(grep -c '>' $R/true.fasta); nbb=$(ls $R/inputs/backbones | wc -l); bsz=$(grep -c '>' $R/inputs/backbones/backbone_1_mafft.txt)
printf '%s\tprotein\t%s\t%s\t%s\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..10}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..10}.fa)\n' \
  $n $nseq $nbb $bsz $(stat -c %s $B/$n.tar.gz) $(sha256sum $B/$n.tar.gz | cut -d' ' -f1) $n $n $n >> $B/MANIFEST.tsv
