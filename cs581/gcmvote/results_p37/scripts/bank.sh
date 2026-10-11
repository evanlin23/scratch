#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# bank.sh K : bank SIMHIGH_RK merge-only inputs + copy result rows into results_p37
set -eu
K=$1; N=SIMHIGH_R$K; C=/home/user/scratch/cs581
R=/opt/work/gcmtrees/reps/$N; V=/opt/work/gcmvote/reps/$N; B=$C/gcmvote/bank; O=$C/gcmvote/results_p37
st=$(mktemp -d); mkdir $st/$N
cp -rL $R/inputs $R/sets $R/aligned $R/true.fasta $R/unaligned.fasta $R/magus.json $st/$N/
cp /opt/data/sim/SIMHIGH/R$K/tree.nwk $st/$N/true_tree.nwk
cp $V/vote/magus/subsets.json $st/$N/subsets.json
mkdir -p $B; tar czf $B/$N.tar.gz -C $st $N; rm -rf $st
[ -f $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
grep -v "^$N	" $B/MANIFEST.tsv > $B/M.tmp || true; mv $B/M.tmp $B/MANIFEST.tsv
nseq=$(grep -c '>' $R/true.fasta); nb=$(ls $R/inputs/backbones | wc -l); bs=$(grep -c '>' $R/inputs/backbones/backbone_1_mafft.txt)
printf '%s\tprotein\t%s\t%s\t%s\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..%s}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..%s}.fa)\n' \
  $N $nseq $nb $bs $(stat -c %s $B/$N.tar.gz) $(sha256sum $B/$N.tar.gz | cut -d' ' -f1) $N $N $nb $N $nb >> $B/MANIFEST.tsv
mkdir -p $O
: > $O/aln.jsonl; for k in 37 38; do cat /opt/work/gcmtrees/reps/SIMHIGH_R$k/results.jsonl /opt/work/gcmvote/reps/SIMHIGH_R$k/results.jsonl >> $O/aln.jsonl 2>/dev/null || true; done
cat /opt/work/gcmtrees/reps/SIMHIGH_R{37,38}/trees.jsonl > $O/trees.jsonl 2>/dev/null || true
grep -h '"SIMHIGH_R3[78]"' /opt/work/gcmtrees/magus.jsonl | grep '"magus"' > $O/magus.jsonl || true
echo BANKED $N
