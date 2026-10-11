#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh NAME : merge-only inputs of a replicate -> cs581/gcmvote/bank/NAME.tar.gz + MANIFEST.tsv row
set -eu
N=$1; lvl=${N%_R*}; r=${N#*_R}
R=/opt/work/gcmtrees/reps/$N; B=/home/user/scratch/cs581/gcmvote/bank; S=/opt/work/bank/$N
rm -rf $S; mkdir -p $S
cp -r $R/inputs $R/sets $R/aligned $R/true.fasta $R/unaligned.fasta $R/magus.json $S/
cp /opt/data/sim/$lvl/R$r/tree.nwk $S/true_tree.nwk
cp /opt/work/gcmvote/reps/$N/vote/magus/subsets.json $S/subsets.json
mkdir -p $B
tar czf $B/$N.tar.gz -C /opt/work/bank $N
M=$B/MANIFEST.tsv
[ -f $M ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $M
grep -v "^$N	" $M > $M.tmp || true; mv $M.tmp $M
nseq=$(grep -c '>' $S/unaligned.fasta); nbb=$(ls $S/inputs/backbones | wc -l)
bsz=$(grep -c '>' $S/sets/s0/backbone_1.fa)
printf '%s\tprotein\t%s\t%s\t%s\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..%s}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..%s}.fa)\n' \
  $N $nseq $nbb $bsz $(stat -c%s $B/$N.tar.gz) $(sha256sum $B/$N.tar.gz | cut -c1-64) $N $N $nbb $N $nbb >> $M
cat $M
