#!/bin/bash
# Backbone-aligner benchmark worker (gcmx.bbtool_bench): MAGUS vs MAGUS with its GCM backbones
# aligned by other tools, paired (merge-only on the same draw) and end to end (timed).
#   bash cs581/code/fanout/bb_worker.sh JOBFILE BRANCH [BASE_BRANCH]
# Phase 1 (clustalo, mafft-auto; draws 0,1,2, draw-major), then phase 2 adds the slow MAFFT
# controls (linsi-noep, ginsi) on the same draws. Run on an otherwise idle machine. Commits and
# pushes cs581/experiments/bbtool/<BRANCH basename>[_p2].jsonl plus per-draw state after every job.
# Before every job BASE_BRANCH (default: the orchestrating session's branch) is merged in, so code
# fixes pushed there take effect.
set -u
if [ -z "${WORKER_REEXEC:-}" ]; then
  # run from a private copy, since the merge below may rewrite this file
  export WORKER_REPO
  WORKER_REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
  copy=$(mktemp /tmp/bb_worker.XXXXXX.sh)
  cp "$0" "$copy"
  WORKER_REEXEC=1 exec bash "$copy" "$(realpath "$1")" "${@:2}"
fi
REPO=$WORKER_REPO
JOBS=$1
BRANCH=$2
BASE=${3:-claude/charming-pasteur-yl6k2v}
DRAWS=${DRAWS:-0,1,2}
OUTDIR=$REPO/cs581/experiments/bbtool
B=$(basename "$BRANCH")
mkdir -p "$OUTDIR/${B}_state" /opt/runs/bbtool
cd "$REPO/cs581/code"
GIT() { git -C "$REPO" -c user.name="Evan Lin" -c user.email="113861384+evanlin23@users.noreply.github.com" "$@"; }
# leftovers of a killed previous invocation (bracket patterns: do not match this shell)
pkill -f "[g]cmx\.(bbtool_bench|run_magus)" 2>/dev/null; pkill -f "[m]afftdir" 2>/dev/null
pkill -x mcl 2>/dev/null; pkill -x clustalo 2>/dev/null; pkill -f "[F]astTree" 2>/dev/null; sleep 2
for phase in 1 2; do
  if [ $phase = 1 ]; then TOOLS=clustalo,mafft-auto; OUT=$OUTDIR/$B.jsonl
  else TOOLS=clustalo,mafft-auto,linsi-noep,ginsi; OUT=$OUTDIR/${B}_p2.jsonl; fi
  for draw in ${DRAWS//,/ }; do
    while read -r name rest; do
      [ -z "${name:-}" ] && continue
      grep -q "\"dataset\": \"$name\", \"draw\": $draw," "$OUT" 2>/dev/null && continue
      GIT pull -q --no-rebase --no-edit origin "$BASE" < /dev/null || echo "$(date +%T) could not merge $BASE"
      echo "$(date +%T) phase $phase draw $draw start $name"
      grep "^$name " "$JOBS" > /tmp/bb_job.txt
      python3 -m gcmx.bbtool_bench /tmp/bb_job.txt "$OUT" /opt/runs/bbtool --draws "$draw" --tools "$TOOLS" \
        --e2e clustalo --threads 4 < /dev/null || echo "$(date +%T) FAILED $name draw $draw"
      cp "/opt/runs/bbtool/${name}_d$draw/state.json" "$OUTDIR/${B}_state/${name}_d$draw.json" 2>/dev/null
      GIT add "$OUTDIR" && GIT commit -q -m "bbtool benchmark: $name draw $draw (phase $phase)"
      for i in 1 2 3 4 5; do GIT push -q -u origin "HEAD:$BRANCH" && break; sleep $((2 ** i)); done
    done < "$JOBS"
  done
done
echo "ALL JOBS DONE"
