#!/bin/bash
# Measured end-to-end benchmark worker: PASTA vs MAGUS vs soft MAGUS (gcmx.e2e_bench).
#   bash cs581/code/fanout/e2e_worker.sh JOBFILE BRANCH
# Run on an otherwise idle machine. Job by job (each restartable per method), commits +
# pushes cs581/experiments/e2e/<BRANCH basename>.jsonl after every job.
set -u
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
JOBS=$(realpath "$1")
BRANCH=$2
OUTDIR=$REPO/cs581/experiments/e2e
OUT=$OUTDIR/$(basename "$BRANCH").jsonl
mkdir -p "$OUTDIR" /opt/runs/e2e
cd "$REPO/cs581/code"
GIT() { git -C "$REPO" -c user.name="Evan Lin" -c user.email="113861384+evanlin23@users.noreply.github.com" "$@"; }
# leftovers of a killed previous invocation (bracket patterns: do not match this shell)
pkill -f "[g]cmx\.(e2e_bench|run_magus|extend|split)" 2>/dev/null; pkill -f "[r]un_pasta\.py" 2>/dev/null
pkill -f "[m]afftdir" 2>/dev/null; pkill -x mcl 2>/dev/null; pkill -x hmmalign 2>/dev/null; pkill -f "[F]astTree" 2>/dev/null
pkill -f "[o]pal\.jar" 2>/dev/null; sleep 2
while read -r name rest; do
  [ -z "${name:-}" ] && continue
  grep -q "\"dataset\": \"$name\"" "$OUT" 2>/dev/null && continue
  echo "$(date +%T) start $name"
  grep "^$name " "$JOBS" > /tmp/e2e_job.txt
  python3 -m gcmx.e2e_bench /tmp/e2e_job.txt "$OUT" /opt/runs/e2e --threads 4 < /dev/null || echo "$(date +%T) FAILED $name"
  GIT add "$OUT" && GIT commit -q -m "e2e benchmark: $name" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
  for i in 1 2 3 4 5; do GIT push -q -u origin "HEAD:$BRANCH" && break; sleep $((2 ** i)); done
done < "$JOBS"
echo "ALL JOBS DONE"
