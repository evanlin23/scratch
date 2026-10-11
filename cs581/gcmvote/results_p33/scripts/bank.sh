#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh NAME : merge-only inputs tarball + MANIFEST row
set -eu
N=$1; lvl=${N%_R*}; r=${N#*_R}
G=/opt/work/gcmtrees/reps/$N; B=/home/user/scratch/cs581/gcmvote/bank; S=/opt/work/bankstage
rm -rf $S/$N; mkdir -p $S/$N $B
cp -r $G/inputs $G/sets $G/aligned $G/true.fasta $G/unaligned.fasta $G/magus.json $S/$N/
cp /opt/data/sim/$lvl/R$r/tree.nwk $S/$N/true_tree.nwk
cp /opt/work/gcmvote/reps/$N/vote/magus/subsets.json $S/$N/subsets.json
tar czf $B/$N.tar.gz -C $S $N
sz=$(stat -c %s $B/$N.tar.gz); sha=$(sha256sum $B/$N.tar.gz | cut -d' ' -f1)
[ -f $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
grep -v "^$N	" $B/MANIFEST.tsv > $B/M.tmp || true; mv $B/M.tmp $B/MANIFEST.tsv
printf '%s\tprotein\t1000\t10\t200\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..10}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..10}.fa)\n' $N $sz $sha $N $N $N >> $B/MANIFEST.tsv
