#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h10b: ROSE 1000M1 R2. Same steps as h10_fresh.sh + h10_vote.sh (branch claude/cs581-gcmvote-h10):
# one fresh MAGUS draw (paper flags) -> rep dir -> linsi/es3/es4/es5 baselines; vote variants (run.py, B=10,
# 4 merges at a time); FastTree GTR+gamma nRF (h10_trees.py, 4 at a time; byte-identical alignments share one
# tree). Rows -> cs581/gcmvote/results_h10b/. Restartable.
set -u
W=/opt/work/h10b; S=/home/user/scratch; V=$S/cs581/gcmvote/code; G=$S/cs581/gcmgen/code; O=$S/cs581/gcmvote/results_h10b
VARIANTS="hard+mask magus es4 hard hard-bb soft soft-bb soft2 soft4"
TREEV="hard soft soft2 soft4 hard-bb soft-bb"
PY=/opt/mm/root/envs/pasta183/bin/python
r=${REP:-R2}; name=1000M1_$r; D=/opt/data/Datasets/ROSE/1000M1/$r
rep=$W/reps/$name; vr=$W/vreps/$name
mkdir -p $W/reps $W/fresh $W/vreps $O
cd $S/cs581/code
if [ ! -f $rep/true.fasta ]; then
  echo "$name $D/rose.aln.true.fasta 25" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads 4 || exit 1
  python3 $S/cs581/protcons/code/pc.py rep $W/fresh/${name}_d0 $rep
fi
mkdir -p $vr
for x in inputs sets true.fasta unaligned.fasta magus.json; do [ -e $vr/$x ] || ln -s $rep/$x $vr/$x; done
# baselines (linsi needs 4 threads -> run alongside the 9 vote merges, 3 at a time)
python3 $G/gg.py run $rep linsi 'linsi#es3' 'linsi#es4' 'linsi#es5' > $W/gg.log 2>&1 &
printf '%s\n' $VARIANTS | xargs -P 3 -I{} python3 $V/run.py $vr {} > $W/vote.log 2>&1
wait
T=$W/trees/$name; mkdir -p $T
args="true=$rep/true.fasta magus=$W/fresh/${name}_d0/magus.fasta es4=$rep/variants/linsi_es_4/out.fasta"
for v in $TREEV; do args="$args vote_$v=$vr/vote/$v/out.fasta"; done
args="$args vote_hard_mask_masked=$vr/vote/hard+mask/out.masked.fasta"
# one tree per distinct alignment: first method with a given md5 is the representative
declare -A rep_of; : > $T/dups.txt
for a in $args; do m=${a%%=*}; f=${a#*=}; [ -s $f ] || continue
  h=$(md5sum < $f | cut -c1-32)
  if [ -n "${rep_of[$h]:-}" ]; then echo "$m ${rep_of[$h]}" >> $T/dups.txt; else rep_of[$h]=$m; ln -sf $f $T/$m.fasta; fi
done
for m in "${rep_of[@]}"; do echo "$m=$T/$m.fasta"; done | \
  xargs -P 4 -I{} nice -n 5 $PY $V/h10_trees.py $name $D/rose.tt $W/trees.jsonl {}
# copy the representative's row for byte-identical alignments
python3 - $W/trees.jsonl $T/dups.txt <<'PY'
import json, sys
rows = [json.loads(l) for l in open(sys.argv[1])]
have = {r["method"]: r for r in rows}
with open(sys.argv[1], "a") as f:
    for l in open(sys.argv[2]):
        m, src = l.split()
        if m not in have and src in have:
            f.write(json.dumps({**have[src], "method": m, "same_alignment_as": src}) + "\n")
PY
{ sed 's/^{/{"src": "gg", /' $rep/results.jsonl; sed 's/^{/{"src": "vote", /' $vr/results.jsonl; } > $O/aln.jsonl
cp $W/magus.jsonl $O/magus.jsonl; cp $W/trees.jsonl $O/trees.jsonl
echo H10B_DONE
