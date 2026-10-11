#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# collect.sh NAME [...] : rebuild results_p43/{aln,trees,magus}.jsonl from the rep dirs; bank NAME's merge-only inputs.
S=/home/user/scratch/cs581/gcmvote; O=$S/results_p43; B=$S/bank; mkdir -p $B
: > $O/aln.jsonl; : > $O/trees.jsonl; : > $O/magus.jsonl
for name in $(ls -d /opt/work/gcmvote/reps/SIMHIGH_R4[34] 2>/dev/null | xargs -n1 basename); do
  G=/opt/work/gcmtrees/reps/$name; V=/opt/work/gcmvote/reps/$name
  cat $G/results.jsonl $V/results.jsonl >> $O/aln.jsonl
  cat $V/trees.jsonl >> $O/trees.jsonl 2>/dev/null
  grep "\"dataset\": \"$name\"" /opt/work/gcmtrees/magus.jsonl >> $O/magus.jsonl
done
for name in "$@"; do
  G=/opt/work/gcmtrees/reps/$name; lvl=${name%_R*}; r=${name#*_R}; tmp=/opt/work/gcmvote/bank_tmp
  rm -rf $tmp; mkdir -p $tmp/$name
  for f in inputs sets aligned true.fasta unaligned.fasta magus.json subsets.json; do
    [ -e $G/$f ] && cp -rL $G/$f $tmp/$name/; done
  cp /opt/data/sim/$lvl/R$r/tree.nwk $tmp/$name/true_tree.nwk
  tar czf $B/$name.tar.gz -C $tmp $name
  bytes=$(stat -c %s $B/$name.tar.gz); sha=$(sha256sum $B/$name.tar.gz | cut -d' ' -f1)
  [ -f $B/MANIFEST.tsv ] || printf 'name\tdtype\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n' > $B/MANIFEST.tsv
  grep -v "^$name	" $B/MANIFEST.tsv > $B/m.tmp; mv $B/m.tmp $B/MANIFEST.tsv
  printf '%s\tprotein\t1000\t10\t200\tyes\t%s\t%s\t%s/true.fasta\t%s/inputs/backbones/backbone_{1..10}_mafft.txt (MAGUS L-INS-i; also %s/aligned/{linsi,fftns2,fftns2-op3}/s0/backbone_{1..10}.fa)\n' \
    $name $bytes $sha $name $name $name >> $B/MANIFEST.tsv
  rm -rf $tmp
done
