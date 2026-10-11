#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank_rep.sh SIMHIGH_Rk : merge-only inputs tarball + MANIFEST row
set -eu
N=$1; r=${N#*_R}; R=/opt/work/gcmtrees/reps/$N; B=/home/user/scratch/cs581/gcmvote/bank
st=/opt/work/bankstage; rm -rf $st/$N; mkdir -p $st/$N $B
cp -r $R/inputs $R/sets $R/aligned $R/true.fasta $R/unaligned.fasta $R/magus.json $st/$N/
cp /opt/data/sim/SIMHIGH/R$r/tree.nwk $st/$N/true_tree.nwk
tar czf $B/$N.tar.gz -C $st $N
[ -f $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
grep -q "^$N	" $B/MANIFEST.tsv || printf '%s\tprotein\t1000\t10\t200\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..10}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..10}.fa)\n' \
  $N $(stat -c %s $B/$N.tar.gz) $(sha256sum $B/$N.tar.gz | cut -c1-64) $N $N $N >> $B/MANIFEST.tsv
