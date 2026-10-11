#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# IQ-TREE (-m LG+G4 --fast -T 2 -seed 1) on true / magus / es4 / hard-bb for SIMHIGH replicates; nRF vs the
# true tree -> cs581/iqtree/results_e/iqtree.jsonl; commit + push after each replicate. Restartable.
#   bash run_e.sh SIMHIGH_R9 SIMHIGH_R10      (reps unpacked under /opt/work/iqe/reps)
set -u
REPO=/home/user/scratch; W=/opt/work/iqe; OUT=$REPO/cs581/iqtree/results_e/iqtree.jsonl
IQ=/opt/mm/root/envs/bio/bin/iqtree3; PY=/opt/mm/root/envs/pasta183/bin/python
VER=$($IQ --version | grep -o 'version [0-9.]*' | head -1 | cut -d' ' -f2)
METHODS="true magus es4 hard-bb"
iq() {  # rep method
  local rep=$1 m=$2 d=$W/iq/$1; mkdir -p $d
  grep -q "\"rep\": \"$rep\", \"method\": \"$m\"" $OUT 2>/dev/null && return
  local aln=$W/reps/$rep/vote/$m/out.fasta; [ $m = true ] && aln=$W/reps/$rep/true.fasta
  if [ ! -s $d/$m.treefile ]; then
    local s=$(date +%s.%N)
    $IQ -s $aln -m LG+G4 --fast -T 2 -seed 1 -pre $d/$m -redo > $d/$m.stdout 2>&1 || { echo "IQFAIL $rep $m"; return; }
    echo "$(date +%s.%N) - $s" | bc > $d/$m.wall
  fi
  local wall=$(cat $d/$m.wall)
  awk -v w=$wall 'BEGIN{if (w > 5400) print "NOTE: over 90 min"}'
  $PY $REPO/cs581/iqtree/code/nrf_e.py $W/reps/$rep/true_tree.nwk $d/$m.treefile $OUT $rep $m $wall $VER
}
for rep in "$@"; do
  echo "== $rep $(date)"
  (cd $REPO/cs581/gcmvote/code && python3 run.py $W/reps/$rep magus es4 hard-bb)
  cat $W/reps/$rep/results.jsonl
  pids=(); for m in $METHODS; do
    iq $rep $m & pids+=($!)
    [ ${#pids[@]} -ge 2 ] && { wait ${pids[0]}; pids=(${pids[@]:1}); }
  done; wait
  mkdir -p $REPO/cs581/iqtree/results_e/$rep
  cp $W/reps/$rep/results.jsonl $REPO/cs581/iqtree/results_e/$rep/vote_results.jsonl
  cp $W/iq/$rep/*.treefile $REPO/cs581/iqtree/results_e/$rep/ 2>/dev/null
  cd $REPO && git add cs581/iqtree && git -c user.name="Evan Lin" -c user.email="113861384+evanlin23@users.noreply.github.com" \
    commit -q -m "iqtree_e: IQ-TREE nRF for $rep (true/magus/es4/hard-bb)

Co-Authored-By: Claude <noreply@anthropic.com>" && \
    for i in 1 2 3 4; do git push -q -u origin claude/cs581-iqtree-e && break; sleep $((2**i)); done
  echo "== done $rep $(date)"
done
echo ALLDONE
