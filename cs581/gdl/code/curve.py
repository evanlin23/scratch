"""Species-tree error vs number of gene families (true gene trees, pure GDL, no ILS)
on 30-taxon Yule trees.  Usage: python curve.py SETTING REP out.jsonl"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gdlsim import yule_tree  # noqa: E402
from probe import run  # noqa: E402

SETTINGS = {
    # name: (lam, mu) for internal branches, (lam, mu) for terminal branches
    "critical": ((2.0, 2.0), (2.0, 2.0)),
    "supercrit": ((2.0, 1.0), (2.0, 1.0)),
    "adversarial": ((3.0, 0.5), (0.0, 3.0)),
}
KS = [int(x) for x in os.environ.get("KS", "25,100,500,2000").split(",")]
setting, rep, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
inner, term = SETTINGS[setting]
tree = yule_tree(30, seed=1000 + rep, height=1.0)


def rule(st, v):
    if st.parent[v] < 0:
        return inner
    return term if not st.children[v] else inner


done = set()
if os.path.exists(out):
    for l in open(out):
        r = json.loads(l)
        done.add((r["setting"], r["rep"], r["nfam"]))
for k in KS:
    if (setting, rep, k) in done:
        continue
    r = run(tree, rule, (0, 0), k, 5000 * rep + k, 4, astral=True)
    r.update({"setting": setting, "rep": rep, "tree": tree})
    with open(out, "a") as f:
        f.write(json.dumps(r) + "\n")
    print(setting, rep, k, {m: v["FN"] for m, v in r["methods"].items()}, flush=True)
