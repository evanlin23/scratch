#!/bin/bash
# Controlled runtime benchmark worker (run on an otherwise idle machine).
#   bash cs581/code/fanout/timing_worker.sh JOBFILE BRANCH
# Runs gcmx.timing_bench job by job; commits + pushes the results file
# cs581/experiments/timing/<BRANCH basename>.jsonl after every job. Restartable.
set -u
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
JOBS=$(realpath "$1")
BRANCH=$2
OUTDIR=$REPO/cs581/experiments/timing
OUT=$OUTDIR/$(basename "$BRANCH").jsonl
mkdir -p "$OUTDIR" /opt/runs/timing
cd "$REPO/cs581/code"
GIT() { git -C "$REPO" -c user.name="Evan Lin" -c user.email="113861384+evanlin23@users.noreply.github.com" "$@"; }
pkill -f "gcmx\.(timing_bench|run_magus|extend)" 2>/dev/null; pkill -f mafftdir 2>/dev/null; sleep 2
while read -r name rest; do
  [ -z "${name:-}" ] && continue
  grep -q "\"dataset\": \"$name\"" "$OUT" 2>/dev/null && continue
  echo "$(date +%T) start $name"
  grep "^$name " "$JOBS" > /tmp/timing_job.txt
  python3 -m gcmx.timing_bench /tmp/timing_job.txt "$OUT" /opt/runs/timing --threads 4 < /dev/null \
    || echo "$(date +%T) FAILED $name"
  GIT add "$OUT" && GIT commit -q -m "timing benchmark: $name" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
  for i in 1 2 3 4 5; do GIT push -q -u origin "HEAD:$BRANCH" && break; sleep $((2 ** i)); done
done < "$JOBS"
echo "ALL JOBS DONE"
