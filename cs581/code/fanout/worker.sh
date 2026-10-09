#!/bin/bash
# Fan-out worker: for each job, run MAGUS with the MAGUS paper's exact settings,
# then the merge-variant pilot on the cached subalignments/backbones, and commit
# + push the (small) results to BRANCH after every job.
#
#   bash cs581/code/fanout/worker.sh JOBFILE BRANCH
#
# JOBFILE lines: NAME SOURCE_TRUE_ALIGNMENT NUM_SUBSETS [PILOT_ONLY]
# (SOURCE is absolute or relative to the repo root.) Restartable: jobs whose
# results are already committed under cs581/experiments/runs/NAME are skipped,
# and leftovers of a killed previous run are cleaned up first.
set -u
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
JOBS=$(realpath "$1")
BRANCH=$2
OUT=$REPO/cs581/experiments/runs
RUNS=/opt/runs
mkdir -p "$OUT" "$RUNS"
cd "$REPO/cs581/code"

# leftovers from a previous (killed) invocation
pkill -f "gcmx\.(prep|pilot|experiment|run_magus)" 2>/dev/null
pkill -f "mafftdir" 2>/dev/null
pkill -x dvtditr 2>/dev/null; pkill -x tbfast 2>/dev/null; pkill -x disttbfast 2>/dev/null
pkill -x mcl 2>/dev/null; pkill -f "FastTree" 2>/dev/null; pkill -x hmmalign 2>/dev/null
sleep 2

push() {
  for i in 1 2 3 4 5; do
    git -C "$REPO" push -q -u origin "HEAD:$BRANCH" && return 0
    sleep $((2 ** i))
  done
  echo "push failed (will retry after next job)"
}

exec 3< "$JOBS"
while read -r -u 3 name src k only; do
  [ -z "${name:-}" ] && continue
  [ -f "$OUT/$name/prep.json" ] && { echo "skip $name (done)"; continue; }
  case "$src" in /*) ;; *) src=$REPO/$src ;; esac
  flags="--maxsubsetsize 0 --maxnumsubsets $k --decompstrategy pastastyle --decompskeletonsize 300
         --graphbuildmethod mafft --graphbuildhmmextend false --graphclustermethod mcl
         --graphtracemethod minclusters --graphtraceoptimize false -r 10 -m 200 -f 4"
  echo "$(date +%T) start $name"
  if [ ! -f "$RUNS/$name/prep.json" ]; then
    python3 -m gcmx.prep "$src" "$RUNS/$name" --threads "$(nproc)" $flags \
      || { echo "$(date +%T) FAILED prep $name"; continue; }
  fi
  PILOT_ONLY=${only:-} PILOT_JOBS=2 python3 -m gcmx.pilot "$RUNS/pilot" "$RUNS/$name" \
    || echo "$(date +%T) pilot errors for $name"

  mkdir -p "$OUT/$name"
  cp "$RUNS/$name/prep.json" "$OUT/$name/"
  python3 -c "
import json, sys
rows = [l for l in open(sys.argv[1]) if json.loads(l)['dataset'] == sys.argv[2]]
open(sys.argv[3], 'w').writelines(rows)" "$RUNS/pilot/results.jsonl" "$name" "$OUT/$name/pilot.jsonl"
  # cached inputs, so later merge experiments need not rerun MAFFT
  tar -C "$RUNS/$name" -cJf "$OUT/$name/inputs.tar.xz" inputs
  git -C "$REPO" add "cs581/experiments/runs/$name"
  git -C "$REPO" -c user.name="Claude" -c user.email="noreply@anthropic.com" commit -q \
    -m "fanout: MAGUS (paper settings) + merge variants on $name" \
    -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
  push
  echo "$(date +%T) done $name"
done
push
echo "ALL JOBS DONE"
