#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# finish.sh K : bank tarball + MANIFEST row + results_p47 rows for SIMHIGH_RK
set -eu
K=$1; N=SIMHIGH_R$K; S=/home/user/scratch/cs581; G=$S/gcmvote
R=/opt/work/gcmtrees/reps/$N; V=/opt/work/gcmvote/reps/$N; O=$G/results_p47; B=$G/bank
mkdir -p $O $B /opt/work/bank/$N
D=/opt/work/bank/$N
rm -rf $D/*; cp -rL $R/inputs $R/sets $R/aligned $R/true.fasta $R/unaligned.fasta $R/magus.json $D/
[ -f $R/subsets.json ] && cp $R/subsets.json $D/
cp /opt/data/sim/SIMHIGH/R$K/tree.nwk $D/true_tree.nwk
tar czf $B/$N.tar.gz -C /opt/work/bank $N
[ -f $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
grep -v "^$N	" $B/MANIFEST.tsv > $B/M.tmp || true; mv $B/M.tmp $B/MANIFEST.tsv
nseq=$(grep -c '>' $R/unaligned.fasta)
printf '%s\tprotein\t%s\t10\t200\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..10}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..10}.fa)\n' \
  $N $nseq $(stat -c %s $B/$N.tar.gz) $(sha256sum $B/$N.tar.gz | cut -d' ' -f1) $N $N $N >> $B/MANIFEST.tsv
# results: per-rep files under /opt/work/p47, concatenated
P=/opt/work/p47; mkdir -p $P
cat $R/results.jsonl $V/results.jsonl > $P/$N.aln.jsonl
cp $R/trees.jsonl $P/$N.trees.jsonl
grep "\"dataset\": \"$N\"" /opt/work/gcmtrees/magus.jsonl > $P/$N.magus.jsonl
for f in aln trees magus; do cat $P/*.$f.jsonl > $O/$f.jsonl; done
echo FINISH_OK $N
