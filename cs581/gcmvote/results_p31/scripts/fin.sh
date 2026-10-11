#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# fin.sh NAME : bank tarball + MANIFEST row, results_p31 jsonl rows, magus-reproduction check
set -eu
n=$1; lvl=${n%_R*}; r=${n#*_R}
S=/home/user/scratch/cs581; O=$S/gcmvote/results_p31; B=$S/gcmvote/bank
R=/opt/work/gcmtrees/reps/$n; V=/opt/work/gcmvote/reps/$n
mkdir -p $B /opt/work/bank/$n
st=/opt/work/bank/$n
rm -rf $st/*; mkdir -p $st/$n
for f in inputs sets aligned true.fasta unaligned.fasta magus.json subsets.json; do
  [ -e $R/$f ] && cp -rL $R/$f $st/$n/; done
cp /opt/data/sim/$lvl/R$r/tree.nwk $st/$n/true_tree.nwk
tar czf $B/$n.tar.gz -C $st $n
sz=$(stat -c %s $B/$n.tar.gz); sha=$(sha256sum $B/$n.tar.gz | cut -d' ' -f1)
nb=$(ls $R/inputs/backbones | wc -l)
[ -f $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
grep -v "^$n	" $B/MANIFEST.tsv > $B/M.tmp || true; mv $B/M.tmp $B/MANIFEST.tsv
printf '%s\tprotein\t1000\t%s\t200\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..10}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..10}.fa)\n' \
  $n $nb $sz $sha $n $n $n >> $B/MANIFEST.tsv
grep -h "\"$n\"" $R/results.jsonl | grep -v '"B"' >> $O/aln.jsonl
cat $V/results.jsonl >> $O/aln.jsonl
cat $R/trees.jsonl >> $O/trees.jsonl
grep -h "\"dataset\": \"$n\"" /opt/work/gcmtrees/magus.jsonl >> $O/magus.jsonl
echo "check magus reproduction:"; python3 $O/scripts/norm.py $R/variants/linsi/out.fasta $V/vote/magus/out.fasta
