#!/bin/bash
# End-to-end BSCAMPP runs with a chosen EPA-ng build. Restartable.
# Usage: run_bscampp.sh <datadir> "<variants>" "<sizes>" "<qtypes>"
#   variants: stock fix bug1only bug2only (binaries /opt/tools/epa/epa-ng-<variant>), or pplacer
#   (BSCAMPP(p): bundled pplacer 1.1.alpha19 + taxit; model from $D/RAxML_info.REF, same tree)
#   env: TAG (suffix for the output dir/label, e.g. _s0), PYTHONHASHSEED (BSCAMPP is otherwise
#   not deterministic: set iteration order changes subtree construction)
D=$1; VARS=$2; SIZES=$3; QT=$4
P=/opt/mm/root/envs/place/bin
for v in $VARS; do for b in $SIZES; do for q in $QT; do
  out=$D/bscampp/${v}_b${b}_$q${TAG}
  [ -s $out/result.jplace ] && continue
  rm -rf $out; mkdir -p $out
  if [ $v = pplacer ]; then pm="--placement-method pplacer"; info=$D/RAxML_info.REF; else pm=""; info=$D/rx.raxml.bestModel; fi
  EPA_BIN=/opt/tools/epa/epa-ng-$v EPA_CALLLOG=$out/epa_calls.tsv /usr/bin/time -v -o $out/time.txt \
    $P/run_bscampp.py $pm -i $info -t $D/rx.raxml.bestTree \
     -a $D/backbone.fa -q $D/query_$q.fa -d $out -o result -b $b -V 5 --threads ${THREADS:-4} --cpus-per-job ${CPJ:-4} > $out/log.txt 2>&1
  w=$(grep "Elapsed (wall" $out/time.txt | awk '{print $NF}'); m=$(grep "Maximum resident" $out/time.txt | awk '{print $NF}')
  echo -e "$v\t$b\t$q${TAG}\t$w\t$m" >> $D/bscampp/times.tsv
done; done; done
