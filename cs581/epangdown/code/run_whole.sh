#!/bin/bash
# EPA-ng on the whole backbone (stock vs fixed), as in BSCAMPP paper Exp. 5. Usage: run_whole.sh <datadir> "<variants>" <qtype>
D=$1; VARS=$2; q=$3
for v in $VARS; do
  out=$D/whole/${v}_$q; [ -s $out/epa_result.jplace ] && continue; mkdir -p $out
  /usr/bin/time -v -o $out/time.txt /opt/tools/epa/epa-ng-$v -t $D/rx.raxml.bestTree -s $D/backbone.fa -q $D/query_$q.fa \
    -m $D/rx.raxml.bestModel -w $out -T 4 --redo > $out/log.txt 2>&1
  echo -e "$v\twhole\t$q\t$(grep 'Elapsed (wall' $out/time.txt | awk '{print $NF}')\t$(grep 'Maximum resident' $out/time.txt | awk '{print $NF}')" >> $D/whole/times.tsv
done
