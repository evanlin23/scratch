#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# h4: after run_base.sh finishes a dataset, run every PREREG.md variant (B=10) with gcmvote's run.py on a
# side replicate dir (one per (symlinks to the rep's inputs/true.fasta, own results.jsonl), collect rows, commit, push.
W=/opt/work/h4; R=/home/user/scratch/cs581; OUT=$R/gcmvote/results_h4; GV=$W/gv/cs581/gcmvote/code
# hard+mask dropped for alignment scoring (its out.fasta equals hard), per the orchestrator; 3 parallel lanes
L1="magus soft soft2 frac0.3"; L2="es4 soft-bb soft4 frac0.4"; L3="hard hard-bb frac0.2 frac0.5"
while read -r name src k; do
  [ -z "$name" ] && continue
  until grep -q "BASE_DONE $name\$" $W/base.log; do
    grep -q ALL_BASE_DONE $W/base.log && break; sleep 60
  done
  grep -q "BASE_DONE $name\$" $W/base.log || { echo "SKIP $name (no base)"; continue; }
  rep=$W/reps/$name; vr=$W/vreps/$name
  mkdir -p $vr; ln -sfn $rep/inputs $vr/inputs; ln -sfn $rep/true.fasta $vr/true.fasta
  for L in "$L1" "$L2" "$L3"; do python3 $GV/run.py $vr $L & done; wait
  mkdir -p $OUT/$name
  cp $rep/results.jsonl $OUT/$name/baselines.jsonl
  cp $vr/results.jsonl $OUT/$name/vote.jsonl
  for t in $vr/vote/*/; do tg=$(basename $t); [ -f $t/model.json ] && cp $t/model.json $OUT/$name/model_$tg.json; done
  grep "\"${name}\"" $W/magus.jsonl > $OUT/$name/magus_run.jsonl
  cd $R/..
  git add cs581/gcmvote/results_h4
  git -c user.name="Evan Lin" -c user.email="113861384+evanlin23@users.noreply.github.com" commit -q \
    -m "gcmvote h4: $name rows (MAGUS draw, linsi es3-5 baselines, vote variants)" \
    -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" \
    -m "Claude-Session: https://claude.ai/code/session_01UKXEnuJh7YdekbkauWYR2B"
  for i in 2 4 8 16 0; do git push -q -u origin claude/cs581-gcmvote-h4 && break; [ $i = 0 ] && break; sleep $i; done
  cd - >/dev/null
  echo "VOTE_DONE $name"
done < "$1"
echo ALL_VOTE_DONE
