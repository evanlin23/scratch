#!/bin/bash
# Rep bank: one tarball per merge-only replicate (no merged outputs), so new GCM variants can be run without MAGUS.
#   bash bank.sh REPNAME [...]   -> bank/REPNAME.tar.gz + a line in bank/MANIFEST.tsv
# Contents: inputs/{subalignments,backbones} (MAGUS's 25 L-INS-i subset alignments and 10 L-INS-i backbones),
# aligned/{linsi,fftns2,fftns2-op3}/s0 (backbone sets as aligned by each tool; linsi = inputs/backbones),
# sets/s0 (unaligned backbone sequence sets), true.fasta, true_tree.nwk, unaligned.fasta, magus.json.
# Use: tar xzf bank/REP.tar.gz -C DIR; then  python3 cs581/gcmgen/code/gg.py run DIR/REP VARIANT ...  (cwd cs581/code)
W=/opt/work/gcmtrees; B=/home/user/scratch/cs581/gcmtrees/bank
M=$B/MANIFEST.tsv
[ -f $M ] || printf "name\tdata_type\tnseq\tn_backbones\tbackbone_size\thas_true_tree\tbytes\tsha256\tpaths\n" > $M
for n in "$@"; do
  rep=$W/reps/$n; lvl=${n%_R*}; r=${n#*_R}; tmp=$(mktemp -d)
  mkdir -p $tmp/$n
  cp -r $rep/inputs $rep/sets $rep/true.fasta $rep/unaligned.fasta $rep/magus.json $tmp/$n/
  mkdir -p $tmp/$n/aligned
  for t in linsi fftns2 fftns2-op3; do [ -d $rep/aligned/$t ] && cp -r $rep/aligned/$t $tmp/$n/aligned/; done
  cp /opt/data/sim/$lvl/R$r/tree.nwk $tmp/$n/true_tree.nwk
  tar czf $B/$n.tar.gz -C $tmp $n
  bytes=$(stat -c %s $B/$n.tar.gz)
  if [ $bytes -gt 52428800 ]; then echo "skip $n ($bytes bytes)"; rm $B/$n.tar.gz; rm -rf $tmp; continue; fi
  nseq=$(grep -c '>' $rep/unaligned.fasta); nbb=$(ls $rep/inputs/backbones | wc -l)
  bbs=$(grep -c '>' $rep/inputs/backbones/$(ls $rep/inputs/backbones | head -1))
  sha=$(sha256sum $B/$n.tar.gz | cut -d' ' -f1)
  paths=$(cd $tmp/$n && find . -maxdepth 2 -mindepth 1 | sed 's|^\./||' | sort | grep -v '/.*/' | tr '\n' ',' | sed 's/,$//')
  sed -i "/^$n\t/d" $M
  printf "%s\tprotein (AliSim LG+G4, %s)\t%s\t%s\t%s\tyes\t%s\t%s\t%s\n" $n $lvl $nseq $nbb $bbs $bytes $sha "$paths" >> $M
  rm -rf $tmp
done
