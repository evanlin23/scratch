#!/bin/bash
# Copy small per-dataset files (results, params, model tree) into results/runs, aggregate, commit, push.
set -e
SIM=${1:-/tmp/claude-0/sims}
M=$(cd "$(dirname "$0")/.." && pwd)
for d in "$SIM"/*_n*_r*/; do
  [ -f "$d/results.json" ] || continue
  n=$(basename "$d"); mkdir -p "$M/results/runs/$n"
  cp "$d"/{results.json,params.json,model.tre} "$M/results/runs/$n/" 2>/dev/null || true
done
python3 "$M/code/aggregate.py" "$SIM" "$M/results" > /dev/null
cd "$M" && git add code results && git commit -qm "${2:-cs581/models: interim simulation results}

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_016r24PFxd7c5QBpCEq4hEba" && git push -q origin claude/cs581-models
