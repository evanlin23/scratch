#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Rep bank entry: everything a merge-only rerun (gg.py run / gcmvote run.py, vote.py) needs, no merged outputs.
# Usage: bank.sh REP_DIR NAME DATATYPE BACKBONE_SIZE HAS_TRUE_TREE
# Tar holds NAME/{inputs/subalignments, inputs/backbones, true.fasta, unaligned.fasta, magus.json, sets/, aligned/}
# (sets/ = bbe.prep output, aligned/linsi = gg.py's view of MAGUS's backbones; both regenerable, kept for speed).
# vote.py writes subsets.json itself at merge time, so there is none to bank.
# Appends a line to cs581/gcmvote/bank/MANIFEST.tsv; skips (and says so) if the tarball is > 50 MB.
set -e
REP=$1; NAME=$2; DT=$3; BBSIZE=$4; TREE=$5
B=/home/user/scratch/cs581/gcmvote/bank
mkdir -p $B
OUT=$B/$NAME.tar.gz
S=$(mktemp -d)
mkdir $S/$NAME
for f in inputs true.fasta unaligned.fasta magus.json sets aligned; do
  [ -e $REP/$f ] && cp -rL $REP/$f $S/$NAME/
done
tar -C $S -cf - $NAME | gzip -9 > $OUT.tmp
rm -rf $S
bytes=$(stat -c %s $OUT.tmp)
if [ $bytes -gt 52428800 ]; then rm -f $OUT.tmp; echo "SKIPPED $NAME: $bytes bytes > 50 MB"; exit 2; fi
mv $OUT.tmp $OUT
nseq=$(grep -c ">" $REP/unaligned.fasta)
nbb=$(ls $REP/inputs/backbones | grep -c "_mafft.txt$")
sha=$(sha256sum $OUT | cut -d" " -f1)
[ -f $B/MANIFEST.tsv ] || printf "name\tdata_type\tnseq\tn_backbones\tbackbone_size\thas_true_tree\ttar_bytes\tsha256\treference_path\tbackbones_path\n" > $B/MANIFEST.tsv
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $NAME $DT $nseq $nbb $BBSIZE $TREE $bytes $sha \
  "$NAME/true.fasta" "$NAME/inputs/backbones/backbone_*_mafft.txt" >> $B/MANIFEST.tsv
echo "BANKED $NAME $bytes $sha"
