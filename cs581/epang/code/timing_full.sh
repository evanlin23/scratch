#!/bin/bash
# Wall clock of each EPA-ng build/flag on the full 9,000-leaf backbone, all 1,000 queries, 4 threads,
# run one at a time on the otherwise idle machine. Usage: timing_full.sh <datadir>
D=$1; P=/opt/mm/root/envs/place/bin; T=/opt/tools/epa
mkdir -p $D/timing
for q in frag full; do
for cfg in "stock|$P/epa-ng|" "simd|$T/epa-ng-simd|" "shift|$T/epa-ng-shift|" "fix|$T/epa-ng-fix2|" \
           "stock_nopremask|$P/epa-ng|--no-pre-mask" "stock_rsoff|$P/epa-ng|--rate-scalers off"; do
  IFS='|' read name bin fl <<< "$cfg"
  out=$D/timing/${name}_$q; [ -s $out/epa_result.jplace ] && continue; mkdir -p $out
  /usr/bin/time -v -o $out/time.txt $bin -t $D/rx.raxml.bestTree -s $D/backbone.fa -q $D/query_$q.fa \
     -m $D/rx.raxml.bestModel -w $out -T 4 --redo $fl > $out/log.txt 2>&1 || rm -f $out/epa_result.jplace
  w=$(grep "Elapsed (wall" $out/time.txt | awk '{print $NF}'); m=$(grep "Maximum resident" $out/time.txt | awk '{print $NF}')
  echo -e "$name\t$q\t$w\t$m\t$(tail -1 $out/log.txt | cut -c1-60)" >> $D/timing/times.tsv
done; done
