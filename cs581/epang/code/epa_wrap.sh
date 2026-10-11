#!/bin/bash
# Wrapper used as BSCAMPP's epa-ng: injects $EPA_EXTRA flags; if $EPA_REEST=1, re-estimates
# GTR+G model and branch lengths on the subtree with RAxML-NG --evaluate before placement.
P=/opt/mm/root/envs/place/bin
args=("$@")
if [ "${EPA_REEST:-0}" = 1 ]; then
  for ((i=0;i<${#args[@]};i++)); do
    case "${args[$i]}" in -m) mi=$((i+1));; -t) ti=$((i+1));; -s) si=$((i+1));; -w) wi=$((i+1));; esac
  done
  w=${args[$wi]}
  $P/raxml-ng --evaluate --msa "${args[$si]}" --tree "${args[$ti]}" --model GTR+G --prefix "$w/reest" \
     --threads 1 --redo --nofiles interim --log ERROR > /dev/null 2>&1
  if [ -s "$w/reest.raxml.bestModel" ]; then
    args[$mi]="$w/reest.raxml.bestModel"; args[$ti]="$w/reest.raxml.bestTree"
  fi
fi
exec ${EPA_BIN:-$P/epa-ng} "${args[@]}" ${EPA_EXTRA}
