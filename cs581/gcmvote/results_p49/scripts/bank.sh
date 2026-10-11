#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh SIMHIGH_Rk : merge-only inputs tarball + MANIFEST row (same layout as helper t20's bank)
set -eu
N=$1; lvl=${N%_R*}; r=${N#*_R}
R=/opt/work/gcmtrees/reps/$N; B=/home/user/scratch/cs581/gcmvote/bank; S=/opt/work/bank_stage
rm -rf $S/$N; mkdir -p $S/$N $B
cp -r $R/inputs $R/sets $R/aligned $R/true.fasta $R/unaligned.fasta $R/magus.json $S/$N/
[ -e $R/subsets.json ] && cp $R/subsets.json $S/$N/
cp /opt/data/sim/$lvl/R$r/tree.nwk $S/$N/true_tree.nwk
tar czf $B/$N.tar.gz -C $S $N
[ -s $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
nseq=$(grep -c '>' $R/unaligned.fasta); nbb=$(ls $R/inputs/backbones | wc -l)
bsz=$(grep -c '>' $R/inputs/backbones/backbone_1_mafft.txt)
grep -v "^$N	" $B/MANIFEST.tsv > $B/MANIFEST.tmp; mv $B/MANIFEST.tmp $B/MANIFEST.tsv
printf '%s\tprotein\t%s\t%s\t%s\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..%s}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..%s}.fa)\n' \
  $N $nseq $nbb $bsz $(stat -c %s $B/$N.tar.gz) $(sha256sum $B/$N.tar.gz | cut -c1-64) $N $N $nbb $N $nbb >> $B/MANIFEST.tsv
tail -1 $B/MANIFEST.tsv
