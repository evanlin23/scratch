#!/bin/bash
# Collect fan-out results from the worker branches into this checkout and
# regenerate the validation report and the merge-variant summary.
#   bash cs581/code/fanout/aggregate.sh
set -u
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
cd "$REPO"
git fetch -q origin '+refs/heads/claude/cs581-worker-*:refs/remotes/origin/claude/cs581-worker-*'
for ref in $(git for-each-ref --format='%(refname:short)' 'refs/remotes/origin/claude/cs581-worker-*'); do
  if git cat-file -e "$ref:cs581/experiments/runs" 2>/dev/null; then
    git checkout "$ref" -- cs581/experiments/runs
    echo "collected $(git ls-tree --name-only "$ref:cs581/experiments/runs" | wc -l) runs from $ref"
  fi
done
cd cs581/code
python3 -m gcmx.validation_report ../experiments/validation/published_scores.jsonl \
  ../experiments/validation/paper_values.json ../experiments/runs > ../experiments/validation/REPORT.md
python3 -m gcmx.summarize ../experiments/runs > ../experiments/variants/SUMMARY.md
echo "wrote cs581/experiments/validation/REPORT.md and cs581/experiments/variants/SUMMARY.md"
