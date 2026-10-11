#!/bin/bash
# End-to-end BSCAMPP(e) runs on one replicate. Restartable.
# Usage: run_bscampp.sh <datadir> "<variants>" "<sizes>" "<qtypes>"
#   variants: stock (bioconda EPA-ng 0.3.8), fix (patched: epa-ng-fix.patch)
D=$1; VARS=$2; SIZES=$3; QT=$4
P=/opt/mm/root/envs/place/bin
for v in $VARS; do for b in $SIZES; do for q in $QT; do
  out=$D/bscampp/${v}_b${b}_$q
  [ -s $out/result.jplace ] && continue
  mkdir -p $out
  # b>=5000: one EPA-ng job at a time (two 5000-leaf jobs exceed 15 GB)
  cpj=2; [ $b -ge 5000 ] && cpj=4
  case $v in stock) bin=$P/epa-ng;; fix) bin=/opt/tools/epa/epa-ng-fix2;; esac
  EPA_BIN=$bin /usr/bin/time -v -o $out/time.txt $P/run_bscampp.py -i $D/rx.raxml.bestModel -t $D/rx.raxml.bestTree \
     -a $D/backbone.fa -q $D/query_$q.fa -d $out -o result.jplace -b $b -V 5 --threads 4 --cpus-per-job $cpj > $out/log.txt 2>&1
  w=$(grep "Elapsed (wall" $out/time.txt | awk '{print $NF}'); m=$(grep "Maximum resident" $out/time.txt | awk '{print $NF}')
  [ -s $out/result.jplace ] && echo -e "$v\t$b\t$q\t$w\t$m\t$cpj" >> $D/bscampp/times.tsv
done; done; done
