#!/bin/bash
# Fast-MAGUS benchmark worker: runs fm_bench.py job by job on an otherwise idle machine and
# commits + pushes its result file (cs581/fastmagus/results/<NAME>.jsonl) after every job.
#
#   bash cs581/fastmagus/code/fm_worker.sh NAME
#
# Jobs come from cs581/fastmagus/code/jobs/<NAME>.txt on origin/claude/cs581-fastmagus,
# re-read before every job (so jobs/variants can be added while the worker runs).
# Restartable: finished (job, variant) pairs are skipped (state in /opt/runs/fm/).
set -u
NAME=$1
BRANCH=claude/cs581-fastmagus
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
OUT=$REPO/cs581/fastmagus/results/$NAME.jsonl
CODE=$REPO/cs581/fastmagus/code
mkdir -p "$(dirname "$OUT")" /opt/runs/fm
GIT() { git -C "$REPO" -c user.name="Claude" -c user.email="noreply@anthropic.com" "$@"; }
# leftovers of a killed previous invocation (bracket patterns: do not match this shell)
pkill -f "[f]m_bench\.py" 2>/dev/null; pkill -f "[g]cmx\.(run_magus|extend|split)" 2>/dev/null
pkill -f "[m]afftdir" 2>/dev/null; pkill -x mcl 2>/dev/null; pkill -x hmmalign 2>/dev/null
pkill -f "[F]astTree" 2>/dev/null; pkill -x clustalo 2>/dev/null; sleep 2

push() {
  GIT add "$OUT" 2>/dev/null && GIT commit -q -m "fastmagus results: $NAME $1" \
    -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || return 0
  for i in 1 2 3 4 5; do
    GIT pull -q --rebase origin "$BRANCH" && GIT push -q -u origin "HEAD:$BRANCH" && return 0
    sleep $((2 ** i))
  done
  echo "push failed (will retry after the next job)"
}

while true; do
  for i in 1 2 3 4; do GIT fetch -q origin "$BRANCH" && break; sleep $((2 ** i)); done
  GIT show "origin/$BRANCH:cs581/fastmagus/code/jobs/$NAME.txt" > /tmp/fm_jobs_$NAME.txt || exit 1
  ran=0
  while read -r job rest; do
    [ -z "${job:-}" ] && continue
    case "$job" in \#*) continue ;; esac
    want=$(echo "$rest" | awk '{print $3}')
    st=/opt/runs/fm/$job/state.json
    if [ -f "$st" ] && python3 -c "import json,sys; r=json.load(open('$st')); sys.exit(0 if all(v in r or 'error' in r.get(v.split('+')[0], {}) for v in '$want'.split(',')) else 1)"; then
      continue
    fi
    fails=$(cat "/opt/runs/fm/$job.fails" 2>/dev/null || echo 0)
    [ "$fails" -ge 2 ] && continue
    echo "$(date +%T) start $job ($want)"
    grep "^$job " /tmp/fm_jobs_$NAME.txt > /tmp/fm_job_$NAME.txt
    python3 "$CODE/fm_bench.py" /tmp/fm_job_$NAME.txt "$OUT" /opt/runs/fm --threads 4 < /dev/null || { echo "$(date +%T) FAILED $job"; echo $((fails + 1)) > "/opt/runs/fm/$job.fails"; }
    push "$job"
    ran=1
    break  # re-read the job list after every job
  done < /tmp/fm_jobs_$NAME.txt
  [ $ran = 0 ] && break
done
push final
echo "ALL JOBS DONE"
