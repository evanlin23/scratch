#!/bin/bash
# Keep building tree jobs for finished soft-MAGUS alignments (+ TRUE) and run them (1 worker),
# until run_soft.sh is gone and nothing is left. METHOD=ft|iqfast
METHOD=${METHOD:-ft}
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=/opt/runs/soft/out; T=/opt/runs/alncrit/trees; J=/opt/runs/alncrit/jobs_soft_$METHOD.tsv
while true; do
  : > $J
  for d in $OUT/*/; do
    name=$(basename $d); ds=${name%_R*}; r=${name##*_}
    tt=/opt/data/published/$ds/$r/true_tree.tre
    for v in default slow slow-soft-m3; do
      [ -f $d/$v.fasta ] && printf "%s/%s\t%s\t%s\t%s\n" $name $v $d/$v.fasta $tt /opt/runs/soft/trees/$name/$v.$METHOD.tre >> $J
    done
    [ -f $d/slow-soft-m3.fasta ] && printf "%s/true\t%s\t%s\t%s\n" $name /opt/data/published/$ds/$r/true_align.txt $tt $T/$ds/$r/true_align.$METHOD.tre >> $J
  done
  python3 $HERE/runtrees.py $J /opt/runs/alncrit/trees_soft.jsonl --workers 1 --method $METHOD
  pgrep -f run_soft.sh > /dev/null || { python3 $HERE/runtrees.py $J /opt/runs/alncrit/trees_soft.jsonl --workers 1 --method $METHOD; break; }
  sleep 120
done
