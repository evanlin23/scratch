#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h10c: ROSE 1000M1 R3. One fresh MAGUS draw (paper flags) -> rep dir -> linsi/es3/es4/es5 baselines
# (as cs581/gcmgen/code/fresh.sh), then the PREREG vote variants via run.py (4 merges at a time) and FastTree
# GTR+gamma nRF trees with h10_trees.py (h10's scorer; 4 at a time; byte-identical alignments share one tree).
# Rows -> cs581/gcmvote/results_h10c/. Restartable.
W=/opt/work/h10c; S=/home/user/scratch; G=$S/cs581/gcmgen/code; V=$S/cs581/gcmvote/code; O=$S/cs581/gcmvote/results_h10c
VARIANTS="hard+mask magus es4 hard hard-bb soft soft-bb soft2 soft4"
TREEV="hard soft soft2 soft4 hard-bb soft-bb"
PY=/opt/mm/root/envs/pasta183/bin/python
r=R3; name=1000M1_$r; rep=$W/reps/$name; vr=$W/vreps/$name; T=$W/trees/$name
mkdir -p $W/reps $W/fresh $vr $T $O
cd $S/cs581/code
if [ ! -f $rep/results.jsonl ] || [ $(grep -c '"variant"' $rep/results.jsonl) -lt 4 ]; then
  echo "$name /opt/data/Datasets/ROSE/1000M1/$r/rose.aln.true.fasta 25" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads 4 || exit 1
  [ -d $rep/inputs ] || python3 $S/cs581/protcons/code/pc.py rep $W/fresh/${name}_d0 $rep
  python3 $G/gg.py run $rep linsi 'linsi#es3' 'linsi#es4' 'linsi#es5'
fi
for x in inputs sets true.fasta unaligned.fasta magus.json; do [ -e $vr/$x ] || ln -s $rep/$x $vr/$x; done
printf '%s\n' $VARIANTS | xargs -P 4 -I{} python3 $V/run.py $vr {}
args="true=$rep/true.fasta magus=$W/fresh/${name}_d0/magus.fasta es4=$rep/variants/linsi_es_4/out.fasta"
for v in $TREEV; do args="$args vote_$v=$vr/vote/$v/out.fasta"; done
args="$args vote_hard_mask_masked=$vr/vote/hard+mask/out.masked.fasta"
declare -A first; : > $W/dups.txt; uniq=""
for a in $args; do m=${a%%=*}; f=${a#*=}; [ -s $f ] || continue
  h=$(md5sum < $f | cut -c1-32)
  if [ -n "${first[$h]}" ]; then echo "$m ${first[$h]}" >> $W/dups.txt; else first[$h]=$m; ln -sf $f $T/$m.fasta; uniq="$uniq $m=$T/$m.fasta"; fi
done
printf '%s\n' $uniq | xargs -P 4 -I{} nice -n 5 $PY $V/h10_trees.py $name /opt/data/Datasets/ROSE/1000M1/$r/rose.tt $W/trees.jsonl {}
# duplicate alignments: copy the tree row of the identical one, marked
python3 - $W/trees.jsonl $W/dups.txt <<'P'
import json, sys
rows = [json.loads(l) for l in open(sys.argv[1])]
have = {r["method"] for r in rows}
by = {r["method"]: r for r in rows}
with open(sys.argv[1], "a") as f:
    for l in open(sys.argv[2]):
        m, src = l.split()
        if m not in have and src in by:
            r = dict(by[src], method=m, same_alignment_as=src); f.write(json.dumps(r) + "\n")
P
{ sed 's/^{/{"src": "gg", /' $rep/results.jsonl; sed 's/^{/{"src": "vote", /' $vr/results.jsonl; } > $O/aln.jsonl
cp $W/magus.jsonl $O/magus.jsonl; cp $W/trees.jsonl $O/trees.jsonl
echo H10C_DONE
