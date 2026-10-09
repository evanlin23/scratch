"""Re-run every confirmed configuration where true-tag ASTRID-Pro was wrong, with and
without root counting (20000 families). Usage: python recheck_root.py confirm.jsonl out.jsonl"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe import run  # noqa: E402

rows = [json.loads(l) for l in open(sys.argv[1])]
with open(sys.argv[2], "w") as f:
    for r in rows:
        last = r["runs"][-1].get("methods", {})
        if last.get("pro-mean-truetags", {}).get("FN", 0) == 0 and last.get("multi-mean", {}).get("FN", 0) == 0:
            continue
        rates = {k: tuple(v) for k, v in r["rates"].items()}
        x = run(r["tree"], rates, (0, 0), 20000, 31, 2, astral=False,
                modes=[("multi", "mean", False), ("pro", "mean", True), ("pro", "mean", "root")])
        out = {"tree": r["tree"], "rates": r["rates"], "confirm_20k": {k: v["FN"] for k, v in last.items()},
               "recheck": {k: v["FN"] for k, v in x["methods"].items()}}
        f.write(json.dumps(out) + "\n"); f.flush()
        print(out["tree"], out["confirm_20k"], out["recheck"], flush=True)
