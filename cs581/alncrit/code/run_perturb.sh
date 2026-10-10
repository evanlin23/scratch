#!/bin/bash
# split (FN-only, longer) perturbations of the TRUE alignment on the soft replicates, scored + FastTree.
HERE=$(cd "$(dirname "$0")" && pwd); O=/opt/runs/perturb; mkdir -p $O
J=$O/jobs.tsv; : > $J
for d in /opt/runs/soft/out/*/; do name=$(basename $d); ds=${name%_R*}; r=${name##*_}
  T=/opt/data/published/$ds/$r/true_align.txt
  for x in "split 0.1" "split 0.2" "shift 0.1"; do set -- $x; f=$O/$name.$1_$2.fa
    [ -f $f ] || python3 $HERE/perturb.py $T $f $1 $2 --seed 1
    grep -q "\"$f\"" $O/scores.jsonl 2>/dev/null || (cd $HERE/../../code && python3 -c "
import json,sys; from gcmx import score; r=score.fastsp(sys.argv[1],sys.argv[2]); r['aln']=sys.argv[2]; r['key']=sys.argv[3]
open('$O/scores.jsonl','a').write(json.dumps(r)+'\n')" $T $f $name/$1_$2)
    printf "%s/%s_%s\t%s\t%s\t%s\n" $name $1 $2 $f /opt/data/published/$ds/$r/true_tree.tre $O/trees/$name.$1_$2.ft.tre >> $J
  done
done
python3 $HERE/runtrees.py $J $O/trees.jsonl --workers 2 --method ft
echo PERTURB_DONE
