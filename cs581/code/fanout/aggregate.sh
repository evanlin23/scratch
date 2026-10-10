#!/bin/bash
# Collect fan-out results from the worker branches into this checkout and
# regenerate the validation report and the merge-variant summary.
#   bash cs581/code/fanout/aggregate.sh
# Worker branches merge the base branch, so a branch can hold (possibly stale)
# copies of other workers' replicates: for every replicate the copy with the
# most pilot rows wins. Files are written with `git show` into the git-ignored
# cs581/experiments/runs_agg/ (the tracked runs/ copies on the base branch are a
# stale snapshot that workers already merged; they are left untouched).
set -u
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
cd "$REPO"
git fetch -q origin '+refs/heads/claude/cs581-worker-*:refs/remotes/origin/claude/cs581-worker-*'
python3 - <<'EOF'
import os, subprocess
git = lambda *a: subprocess.run(["git"] + list(a), capture_output=True)
refs = git("for-each-ref", "--format=%(refname:short)", "refs/remotes/origin/claude/cs581-worker-*").stdout.decode().split()
best = {}
for ref in refs:
    for path in git("ls-tree", "--name-only", ref + ":cs581/experiments/runs").stdout.decode().split():
        rows = git("show", "{}:cs581/experiments/runs/{}/pilot.jsonl".format(ref, path)).stdout.count(b"\n")
        if path not in best or rows > best[path][1]:
            best[path] = (ref, rows)
for name, (ref, rows) in sorted(best.items()):
    out = os.path.join("cs581/experiments/runs_agg", name)
    os.makedirs(out, exist_ok=True)
    for f in ("prep.json", "pilot.jsonl", "inputs.tar.xz"):
        blob = git("show", "{}:cs581/experiments/runs/{}/{}".format(ref, name, f))
        if blob.returncode == 0:
            open(os.path.join(out, f), "wb").write(blob.stdout)
print("collected {} replicates from {} worker branches".format(len(best), len(refs)))
EOF
cd cs581/code
python3 -m gcmx.validation_report ../experiments/validation/published_scores.jsonl \
  ../experiments/validation/paper_values.json ../experiments/runs_agg > ../experiments/validation/REPORT.md
python3 -m gcmx.summarize ../experiments/runs_agg > ../experiments/variants/SUMMARY.md
echo "wrote cs581/experiments/validation/REPORT.md and cs581/experiments/variants/SUMMARY.md"
