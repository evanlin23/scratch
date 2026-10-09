"""Re-run flagged search configurations with more families (and ASTRAL-Pro) to see
whether a wrong limiting tree persists.  Usage: python confirm.py out.jsonl in1.jsonl ..."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe import run  # noqa: E402

out = sys.argv[1]
KEY = ["multi-mean", "pro-mean", "pro-mean-truetags", "pro-min"]
seen = set()
if os.path.exists(out):
    for l in open(out):
        seen.add(json.loads(l)["tree"] + json.dumps(json.loads(l)["rates"], sort_keys=True))
cands = []
for fn in sys.argv[2:]:
    for l in open(fn):
        r = json.loads(l)
        if any(r["methods"].get(k, {}).get("FN", 0) > 0 for k in KEY):
            cands.append(r)
for r in cands:
    key = r["tree"] + json.dumps(r["rates"], sort_keys=True)
    if key in seen:
        continue
    seen.add(key)
    rates = {k: tuple(v) for k, v in r["rates"].items()}
    rows = []
    for nfam, seed in [(5000, 7), (20000, 8)]:
        try:
            x = run(r["tree"], rates, (0, 0), nfam, seed, 2, astral=True)
        except Exception as e:  # noqa
            x = {"error": str(e)}
        rows.append(x)
    with open(out, "a") as f:
        f.write(json.dumps({"tree": r["tree"], "rates": r["rates"], "shape": r["shape"],
                            "search_FN": {k: v["FN"] for k, v in r["methods"].items()},
                            "runs": rows}) + "\n")
    print(r["tree"], [{k: v["FN"] for k, v in x.get("methods", {}).items()} for x in rows], flush=True)
