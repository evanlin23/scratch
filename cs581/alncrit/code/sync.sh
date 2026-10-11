#!/bin/bash
# copy small result files into the repo and push (run periodically)
R=$(git -C "$(dirname "$0")" rev-parse --show-toplevel); D=$R/cs581/alncrit/results
cp /opt/runs/alncrit/trees_pub.jsonl /opt/runs/alncrit/trees_soft.jsonl /opt/runs/alncrit/alnstats.jsonl $D/ 2>/dev/null
cp /opt/runs/perturb/scores.jsonl $D/perturb_scores.jsonl 2>/dev/null; cp /opt/runs/perturb/trees.jsonl $D/perturb_trees.jsonl 2>/dev/null
cp /opt/runs/soft/out/results.jsonl $D/soft_alignments.jsonl; cp /opt/runs/soft/out/timing.jsonl $D/soft_timing.jsonl
cd $R && git add cs581/alncrit/results && git -c user.name=Claude -c user.email=noreply@anthropic.com commit -q -m "alncrit: interim results $(date +%H:%M)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01NgvCDj8KAVSiKtmcxkoUwQ" && for i in 1 2 3 4; do git push -q origin claude/cs581-alncrit && break; sleep $((2**i)); done
