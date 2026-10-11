#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh NAME : tar the merge-only inputs of the gcmtrees rep (layout of banks t18-t20) and append a MANIFEST row
set -eu
N=$1; lvl=${N%_R*}; r=${N#*_R}
R=/opt/work/gcmtrees/reps/$N; B=/home/user/scratch/cs581/gcmvote/bank; S=/opt/work/bank_stage
rm -rf $S/$N; mkdir -p $S/$N $B
cp -rL $R/inputs $R/sets $R/aligned $R/true.fasta $R/unaligned.fasta $R/magus.json $S/$N/
cp /opt/data/sim/$lvl/R$r/tree.nwk $S/$N/true_tree.nwk
tar czf $B/$N.tar.gz -C $S $N
[ -s $B/MANIFEST.tsv ] || printf 'name\tdata_type\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
nseq=$(grep -c '>' $R/true.fasta); nbb=$(ls $R/inputs/backbones | wc -l)
bs=$(grep -c '>' $R/inputs/backbones/backbone_1_mafft.txt)
printf '%s\tprotein\t%s\t%s\t%s\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..%s}_mafft.txt (L-INS-i; also aligned/{fftns2,fftns2-op3}/s0/)\n' \
  $N $nseq $nbb $bs $(stat -c %s $B/$N.tar.gz) $(sha256sum $B/$N.tar.gz | cut -d' ' -f1) $N $N $nbb >> $B/MANIFEST.tsv
tail -1 $B/MANIFEST.tsv
