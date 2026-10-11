#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h3b: fresh MAGUS draw (paper flags, as gcmgen/code/fresh.sh) + rep dir, gg.py baselines, vote.py variants
# (run.py scoring, one process per variant, 4 lanes), rows copied to results_h3b/. Restartable.
#   bash pipeline.sh NAME   (NAME = BBA0081 | BBA0117)
set -u
S=/home/user/scratch/cs581; W=/opt/work/gcmvote; name=$1; rep=$W/reps/$name
mkdir -p $W
if [ ! -d $rep/sets/s0 ]; then
  echo "$name $S/data/balibase_clean/RV100_$name.fasta 25" > $W/job_$name.txt
  (cd $S/code && python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/draws --draws 0 --tools '' --e2e '' --union '' --threads 4) || exit 1
  python3 $S/protcons/code/pc.py rep $W/draws/${name}_d0 $rep || exit 1
fi
# baselines (gg.py results.jsonl) in one lane; vote variants in their own rep copies' result files
(python3 $S/gcmgen/code/gg.py run $rep linsi 'linsi#es3' 'linsi#es4' 'linsi#es5' > $W/gg_$name.log 2>&1) &
for v in magus es4 hard soft soft2 soft4 hard-bb soft-bb; do echo $v; done | \
  xargs -P 3 -I{} python3 $S/gcmvote/results_h3b/vote_one.py $rep {}
wait
out=$S/gcmvote/results_h3b; mkdir -p $out
cp $rep/results.jsonl $out/$name.gg_results.jsonl
cat $rep/vote_rows/*.json > $out/$name.results.jsonl
for v in magus es4 hard soft soft2 soft4 hard-bb soft-bb; do cp $rep/vote/$v/model.json $out/$name.$v.model.json; done
echo PIPE_DONE $name
