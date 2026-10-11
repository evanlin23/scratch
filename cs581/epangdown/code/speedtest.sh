#!/bin/bash
# Standalone EPA-ng wall-clock on one k-leaf subtree, one thread, each build in turn (sequential).
# Usage: speedtest.sh <datadir> <k> <nq> [reps=1]
D=$1; K=$2; NQ=$3; R=${4:-1}
P=/opt/mm/root/envs/place/bin; C=$(cd "$(dirname "$0")" && pwd)
S=$D/speed/k$K; [ -s $S/tree.nwk ] || $P/python $C/make_subtree.py $D $K $S $NQ
for rep in $(seq 1 $R); do
for q in frag full; do
for v in stock bug1only bug2only fix stock_rsoff; do
  bin=/opt/tools/epa/epa-ng-${v%_rsoff}; extra=""; [ "${v%_rsoff}" != "$v" ] && extra="--rate-scalers off"
  out=$S/${v}_${q}_r$rep; mkdir -p $out
  /usr/bin/time -v -o $out/time.txt $bin -t $S/tree.nwk -s $S/ref.fa -q $S/q_$q.fa -m $D/rx.raxml.bestModel \
     -w $out -T 1 --redo $extra > $out/log.txt 2>&1
  echo -e "$K\t$v\t$q\t$rep\t$(grep 'Elapsed (wall' $out/time.txt | awk '{print $NF}')\t$(grep 'Maximum resident' $out/time.txt | awk '{print $NF}')" >> $D/speed/times.tsv
done; done; done
