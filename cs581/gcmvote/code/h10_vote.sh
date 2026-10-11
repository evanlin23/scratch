#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h10: after h10_fresh.sh finishes a replicate's gg baselines, run every PREREG vote variant (vote.py via
# run.py, B=10; 4 merges at a time) and FastTree GTR+gamma nRF trees (4 at a time); copy rows to
# cs581/gcmvote/results_h10/. Restartable. vote.py/run.py come from branch claude/cs581-gcmvote, checked out as a
# worktree at /opt/gcmvote_wt. hard+mask is run only for its masked tree alignment (its unmasked output = hard).
W=/opt/work/h10; V=/opt/gcmvote_wt/cs581/gcmvote/code; S=/home/user/scratch; O=$S/cs581/gcmvote/results_h10
VARIANTS="hard+mask magus es4 hard hard-bb soft soft-bb soft2 soft4 frac0.2 frac0.3 frac0.4 frac0.5"
TREEV="hard soft soft2 soft4 hard-bb soft-bb"
PY=/opt/mm/root/envs/pasta183/bin/python
mkdir -p $O $W/vreps
for r in ${REPS:-R0 R1}; do
  name=1000M1_$r; rep=$W/reps/$name; vr=$W/vreps/$name
  until [ -f $rep/results.jsonl ] && [ $(grep -c '"variant"' $rep/results.jsonl) -ge 4 ]; do sleep 60; done
  mkdir -p $vr
  for x in inputs sets true.fasta unaligned.fasta magus.json; do [ -e $vr/$x ] || ln -s $rep/$x $vr/$x; done
  printf '%s\n' $VARIANTS | xargs -P 4 -I{} python3 $V/run.py $vr {}
  T=$W/trees/$name; mkdir -p $T
  args="true=$rep/true.fasta magus=$W/fresh/${name}_d0/magus.fasta es4=$rep/variants/linsi_es_4/out.fasta"
  for v in $TREEV; do args="$args vote_$v=$vr/vote/$v/out.fasta"; done
  args="$args vote_hard_mask_masked=$vr/vote/hard+mask/out.masked.fasta"
  for a in $args; do m=${a%%=*}; f=${a#*=}; [ -s $f ] && ln -sf $f $T/$m.fasta; done
  for a in $args; do m=${a%%=*}; [ -e $T/$m.fasta ] && echo "$m=$T/$m.fasta"; done | \
    xargs -P 4 -I{} nice -n 5 $PY $S/cs581/gcmvote/code/h10_trees.py $name /opt/data/Datasets/ROSE/1000M1/$r/rose.tt $W/trees.jsonl {}
  { cat $W/reps/*/results.jsonl | sed 's/^{/{"src": "gg", /'; cat $W/vreps/*/results.jsonl | sed 's/^{/{"src": "vote", /'; } > $O/aln.jsonl
  cp $W/magus.jsonl $O/magus.jsonl; cp $W/trees.jsonl $O/trees.jsonl
  cd $S && git add $O && git -c user.name="Evan Lin" -c user.email="113861384+evanlin23@users.noreply.github.com" \
    commit -q -m "gcmvote h10: 1000M1 $r rows" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_0183UYDvTv7komXdMwUC7GQ6" && \
    for i in 1 2 3 4; do git push -q origin claude/cs581-gcmvote-h10 && break; sleep $((2**i)); done
done
echo H10_VOTE_DONE
