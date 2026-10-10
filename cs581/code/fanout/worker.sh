#!/bin/bash
# Fan-out worker: for each job, run MAGUS with the MAGUS paper's exact settings,
# then the merge-variant pilot on the cached subalignments/backbones, and commit
# + push the (small) results to BRANCH after every job.
#
#   bash cs581/code/fanout/worker.sh JOBFILE BRANCH [BASE_BRANCH]
#
# JOBFILE lines: NAME SOURCE_TRUE_ALIGNMENT NUM_SUBSETS [PILOT_ONLY]
# (SOURCE is absolute or relative to the repo root.)
# - Restartable: jobs with results committed under cs581/experiments/runs/NAME
#   are skipped, and leftovers of a killed previous invocation are cleaned up.
# - Before every job BASE_BRANCH (default: the orchestrating session's branch)
#   is merged in, so code fixes and job-list changes pushed there take effect.
set -u
if [ -z "${WORKER_REEXEC:-}" ]; then
  # run from a private copy, since the merge below may rewrite this file
  export WORKER_REPO
  WORKER_REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
  copy=$(mktemp /tmp/worker.XXXXXX.sh)
  cp "$0" "$copy"
  WORKER_REEXEC=1 exec bash "$copy" "$(realpath "$1")" "${@:2}"
fi
REPO=$WORKER_REPO
JOBS=$1
BRANCH=$2
BASE=${3:-claude/charming-pasteur-yl6k2v}
OUT=$REPO/cs581/experiments/runs
RUNS=/opt/runs
mkdir -p "$OUT" "$RUNS"
cd "$REPO/cs581/code"
GIT() { git -C "$REPO" -c user.name="Evan Lin" -c user.email="113861384+evanlin23@users.noreply.github.com" "$@"; }

# leftovers from a previous (killed) invocation
pkill -f "gcmx\.(prep|pilot|experiment|run_magus)" 2>/dev/null
pkill -f "mafftdir" 2>/dev/null
pkill -x dvtditr 2>/dev/null; pkill -x tbfast 2>/dev/null; pkill -x disttbfast 2>/dev/null
pkill -x mcl 2>/dev/null; pkill -f "FastTree" 2>/dev/null; pkill -x hmmalign 2>/dev/null
sleep 2

push() {
  for i in 1 2 3 4 5; do
    GIT push -q -u origin "HEAD:$BRANCH" && return 0
    sleep $((2 ** i))
  done
  echo "push failed (will retry after next job)"
}

next_job() {
  while read -r name rest; do
    [ -z "${name:-}" ] && continue
    [ -f "$OUT/$name/prep.json" ] || [ -f "$RUNS/FAILED_$name" ] || { echo "$name $rest"; return; }
  done < "$JOBS"
}

while true; do
  GIT pull -q --no-rebase --no-edit origin "$BASE" || echo "$(date +%T) could not merge $BASE"
  job=$(next_job)
  [ -z "$job" ] && break
  read -r name src k only <<< "$job"
  case "$src" in /*) ;; *) src=$REPO/$src ;; esac
  flags="--maxsubsetsize 0 --maxnumsubsets $k --decompstrategy pastastyle --decompskeletonsize 300
         --graphbuildmethod mafft --graphbuildhmmextend false --graphclustermethod mcl
         --graphtracemethod minclusters --graphtraceoptimize false -r 10 -m 200 -f 4"
  echo "$(date +%T) start $name"
  if [ ! -f "$RUNS/$name/prep.json" ]; then
    if ! python3 -m gcmx.prep "$src" "$RUNS/$name" --threads "$(nproc)" $flags < /dev/null; then
      echo "$(date +%T) FAILED prep $name"
      touch "$RUNS/FAILED_$name"
      continue
    fi
  fi
  PILOT_ONLY=${only:-} PILOT_JOBS=2 python3 -m gcmx.pilot "$RUNS/pilot" "$RUNS/$name" < /dev/null \
    || echo "$(date +%T) pilot errors for $name"

  mkdir -p "$OUT/$name"
  cp "$RUNS/$name/prep.json" "$OUT/$name/"
  python3 -c "
import json, sys
rows = [l for l in open(sys.argv[1]) if json.loads(l)['dataset'] == sys.argv[2]]
open(sys.argv[3], 'w').writelines(rows)" "$RUNS/pilot/results.jsonl" "$name" "$OUT/$name/pilot.jsonl"
  # cached inputs, so later merge experiments need not rerun MAFFT
  tar -C "$RUNS/$name" -cJf "$OUT/$name/inputs.tar.xz" inputs
  GIT add "cs581/experiments/runs/$name"
  GIT commit -q -m "fanout: MAGUS (paper settings) + merge variants on $name" \
    -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
  push
  echo "$(date +%T) done $name"
done

# Backfill: rerun the pilot on finished replicates this machine computed, so they
# also get variants added after they were first processed (done variants are skipped).
while read -r name rest; do
  [ -f "$RUNS/$name/prep.json" ] && [ -f "$OUT/$name/prep.json" ] || continue
  GIT pull -q --no-rebase --no-edit origin "$BASE" || true
  only=$(echo "$rest" | awk '{print $3}')
  PILOT_ONLY=${only:-} PILOT_JOBS=2 python3 -m gcmx.pilot "$RUNS/pilot" "$RUNS/$name" < /dev/null \
    || echo "$(date +%T) backfill pilot errors for $name"
  python3 -c "
import json, sys
rows = [l for l in open(sys.argv[1]) if json.loads(l)['dataset'] == sys.argv[2]]
open(sys.argv[3], 'w').writelines(rows)" "$RUNS/pilot/results.jsonl" "$name" "$OUT/$name/pilot.jsonl"
  if ! git -C "$REPO" diff --quiet -- "cs581/experiments/runs/$name"; then
    GIT add "cs581/experiments/runs/$name"
    GIT commit -q -m "fanout: backfill merge variants on $name" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
    push
    echo "$(date +%T) backfilled $name"
  fi
done < "$JOBS"
push
echo "ALL JOBS DONE"
