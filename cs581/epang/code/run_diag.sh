#!/bin/bash
# Diagnostic: CPU (non-SIMD) kernels without rate scalers, via the EPA_FORCE_CPU switch
# of the diag build (epa-ng-diag.patch). Mode dir: diagcpu_<qtype>. Usage: run_diag.sh <datadir> <qtype> "<sizes>"
D=$1; q=$2; SIZES=${3:-"500 1000 2000"}
for k in $SIZES; do for kd in $D/nested/c*/k$k; do
  out=$kd/diagcpu_$q; [ -s $out/epa_result.jplace ] && continue; mkdir -p $out
  s=$(date +%s.%N)
  EPA_FORCE_CPU=1 /opt/tools/epa/epa-ng-diag -t $kd/tree.nwk -s $kd/ref.fa -q $kd/$q.fa -m $D/rx.raxml.bestModel \
     -w $out -T 1 --redo --rate-scalers off > $out/log.txt 2>&1 || rm -f $out/epa_result.jplace
  echo -e "$kd\tdiagcpu\t$q\t$(echo "$(date +%s.%N) - $s" | bc)\t0" >> $D/nested/times.tsv
done; done
