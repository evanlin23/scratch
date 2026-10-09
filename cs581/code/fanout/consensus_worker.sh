#!/bin/bash
# Fan-out worker for the consensus-MSA significance test on the MAGUS paper's published alignments.
#
#   bash cs581/code/fanout/consensus_worker.sh REPLIST BRANCH [BASE_BRANCH]
#
# REPLIST lines: DATASET/REP as in the paper's Results.zip (e.g. 1000M2/R0, balibase/RV100_BBA0039).
# Runs gcmx.consensus_eval in chunks of 6 replicates, commits + pushes the results
# file cs581/experiments/consensus/<BRANCH basename>.jsonl after each chunk.
# Restartable: finished (replicate, scenario) pairs are skipped.
set -u
if [ -z "${WORKER_REEXEC:-}" ]; then
  export WORKER_REPO
  WORKER_REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
  copy=$(mktemp /tmp/cworker.XXXXXX.sh)
  cp "$0" "$copy"
  WORKER_REEXEC=1 exec bash "$copy" "$(realpath "$1")" "${@:2}"
fi
REPO=$WORKER_REPO
LIST=$1
BRANCH=$2
BASE=${3:-claude/charming-pasteur-yl6k2v}
OUTDIR=$REPO/cs581/experiments/consensus
OUT=$OUTDIR/$(basename "$BRANCH").jsonl
WORK=/opt/runs/consensus
mkdir -p "$OUTDIR" "$WORK"
cd "$REPO/cs581/code"
GIT() { git -C "$REPO" -c user.name="Claude" -c user.email="noreply@anthropic.com" "$@"; }

pkill -f "gcmx\.(consensus_eval|run_magus)" 2>/dev/null
pkill -x mcl 2>/dev/null
sleep 2

push() {
  for i in 1 2 3 4 5; do
    GIT push -q -u origin "HEAD:$BRANCH" && return 0
    sleep $((2 ** i))
  done
  echo "push failed (will retry later)"
}

mapfile -t reps < <(grep -v '^\s*$' "$LIST")
for ((i = 0; i < ${#reps[@]}; i += 6)); do
  GIT pull -q --no-rebase --no-edit origin "$BASE" || echo "$(date +%T) could not merge $BASE"
  chunk=("${reps[@]:i:6}")
  echo "$(date +%T) chunk ${chunk[*]}"
  python3 -m gcmx.consensus_eval "$OUT" "$WORK" "${chunk[@]}" --jobs "${CONSENSUS_JOBS:-3}" < /dev/null \
    || echo "$(date +%T) errors in chunk ${chunk[*]}"
  rm -rf "${WORK:?}"/*
  if [ -f "$OUT" ]; then
    GIT add "$OUT"
    GIT commit -q -m "consensus test: ${chunk[*]}" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" && push
  fi
done
push
echo "ALL JOBS DONE"
