#!/bin/bash
# MAGUS (paper settings, merge-only rerun on cached inputs = "default"), MAGUS(Slow) ("slow")
# and slow-soft-m3 on replicates whose MAGUS inputs were cached by the fanout workers.
#   bash run_soft.sh REPLIST   (lines: NAME WORKER_BRANCH, e.g. "1000L1_R1 claude/cs581-worker-2")
# Restartable: replicates with all three variants in results.jsonl are skipped.
set -u
LIST=$(realpath "$1")
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
ROOT=/opt/runs/soft; OUT=$ROOT/out; mkdir -p "$ROOT/reps" "$OUT"
cd "$REPO/cs581/code"
while read -r name branch; do
  [ -z "${name:-}" ] && continue
  n=$(python3 -c "
import json,sys
try: rows=[json.loads(l) for l in open('$OUT/results.jsonl')]
except FileNotFoundError: rows=[]
print(len({r['variant'] for r in rows if r.get('dataset')=='$name' and 'avgErr' in r}))")
  [ "$n" -ge 3 ] && continue
  rep=$ROOT/reps/$name; ds=${name%_R*}; r=${name##*_}
  if [ ! -f "$rep/unaligned.fasta" ]; then
    mkdir -p "$rep"
    git -C "$REPO" show "origin/$branch:cs581/experiments/runs/$name/inputs.tar.xz" | tar -C "$rep" -xJ
    git -C "$REPO" show "origin/$branch:cs581/experiments/runs/$name/prep.json" > "$rep/prep.json"
    python3 -c "
import sys; sys.path.insert(0,'.')
from gcmx import fasta
t=fasta.upper(fasta.read(sys.argv[1])); fasta.write(t,sys.argv[2]); fasta.write(fasta.ungap(t),sys.argv[3])
" "/opt/data/published/$ds/$r/true_align.txt" "$rep/true.fasta" "$rep/unaligned.fasta"
  fi
  s=$(date +%s)
  [ -d "$rep/split_m3" ] || python3 -m gcmx.split "$rep/inputs/subalignments" "$rep/split_m3" 3 --method linkage
  t1=$(date +%s)
  if [ ! -f "$rep/ext_backbones.DONE" ]; then
    python3 -m gcmx.extend "$rep/inputs/backbones" "$rep/unaligned.fasta" "$rep/ext_backbones" --jobs 4 && touch "$rep/ext_backbones.DONE"
  fi
  t2=$(date +%s)
  echo "{\"dataset\": \"$name\", \"variant\": \"_timing\", \"split_seconds\": $((t1-s)), \"extend_seconds\": $((t2-t1))}" >> "$OUT/timing.jsonl"
  python3 -m gcmx.experiment --dataset "$name" --true "$rep/true.fasta" \
    --subalignments "$rep/inputs/subalignments" --backbones "$rep/inputs/backbones" --outdir "$OUT" --jobs 3 \
    --variant "default:" --variant "slow:-b $rep/ext_backbones" \
    --variant "slow-soft-m3:-s $rep/split_m3 -b $rep/ext_backbones"
  echo "$(date +%T) done $name"
done < "$LIST"
