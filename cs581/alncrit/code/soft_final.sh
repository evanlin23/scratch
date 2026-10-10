#!/bin/bash
# After all soft-MAGUS alignments exist: FastTree, then IQ-TREE --fast, on MAGUS, slow-soft-m3,
# MAGUS(Slow) and TRUE alignments of each replicate (2 workers). Restartable.
HERE=$(cd "$(dirname "$0")" && pwd)
mk(){ M=$1; J=/opt/runs/alncrit/jobs_soft_$M.tsv; : > $J
  for d in /opt/runs/soft/out/*/; do name=$(basename $d); ds=${name%_R*}; r=${name##*_}
    tt=/opt/data/published/$ds/$r/true_tree.tre
    for v in default slow-soft-m3 slow; do
      [ -f $d/$v.fasta ] && printf "%s/%s\t%s\t%s\t%s\n" $name $v $d/$v.fasta $tt /opt/runs/soft/trees/$name/$v.$M.tre >> $J
    done
    printf "%s/true\t%s\t%s\t%s\n" $name /opt/data/published/$ds/$r/true_align.txt $tt /opt/runs/soft/trees/$name/true.$M.tre >> $J
  done; }
mk ft; python3 $HERE/runtrees.py /opt/runs/alncrit/jobs_soft_ft.tsv /opt/runs/alncrit/trees_soft.jsonl --workers 2 --method ft
mk iqfast; python3 $HERE/runtrees.py /opt/runs/alncrit/jobs_soft_iqfast.tsv /opt/runs/alncrit/trees_soft.jsonl --workers 2 --method iqfast
echo SOFT_ALL_DONE
