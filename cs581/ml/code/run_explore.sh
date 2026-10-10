#!/bin/bash
# After the held-out test queue: (1) redo the UPP arm with the corrected UPP configuration
# (the first UPP alignments put fragments in UPP's backbone; those rows are moved to
# results/test_upp_misconfigured.jsonl and not reported); (2) alternative 1 starting-tree study.
cd "$(dirname "$0")/.."
while pgrep -x xargs > /dev/null; do sleep 60; done
python3 - <<'PY'
import json
rows = [json.loads(l) for l in open("results/test.jsonl")]
bad = [r for r in rows if r["aln"] == "upp"]
with open("results/test_upp_misconfigured.jsonl", "a") as f:
    for r in bad: f.write(json.dumps(r) + "\n")
with open("results/test.jsonl", "w") as f:
    for r in rows:
        if r["aln"] != "upp": f.write(json.dumps(r) + "\n")
print("moved", len(bad), "UPP rows")
PY
rm -f /opt/data/mlcache/*/R*/upp.fasta /opt/data/mlcache/*/R*/upp.clean.fasta
xargs -P 3 -I CMD bash -c CMD < code/explore.cmds
