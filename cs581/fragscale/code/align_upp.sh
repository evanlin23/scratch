#!/bin/bash
# UPP (SEPP 4.5.6, default settings, PASTA backbone) on one replicate's unaligned sequences.
# -M 0.75: "full length" = within 25% of the 3rd quartile of lengths (the default median of all
# sequences lies between fragments and full-length sequences on HF data). Writes
# $MLDATA/DS/R<rep>/upp.fasta (= UPP's masked alignment: insertion columns removed) and upp.time.
#   bash align_upp.sh DS REP [CPUS]
DS=$1; R=$2; X=${3:-1}
D=${MLDATA:-/opt/data/fscache}/$DS/R$R
[ -s $D/upp.fasta ] && exit 0
E=/opt/mm/root/envs/aln/bin
[ -s $D/unaligned.fasta ] || sed "/^>/!s/[-.]//g" $D/true_align.fasta > $D/unaligned.fasta
W=$D/upp_work; rm -rf $W; mkdir -p $W
PATH=$E:$PATH /usr/bin/time -v -o $D/upp.time $E/run_upp.py -s $D/unaligned.fasta -M 0.75 -x $X -d $W -o upp \
  -p $W/tmp > $W/log.txt 2>&1 && cp $W/upp_alignment_masked.fasta $D/upp.fasta
ls -la $D/upp.fasta
