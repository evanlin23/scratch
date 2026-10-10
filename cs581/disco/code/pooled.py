"""Pooled paired comparisons across all estimated-gene-tree experiments (one row per replicate x experiment)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aggregate import load, compare
EXPS = [("DISCO data 101 sp, 50 genes", "/opt/runs/disco", lambda r: "k50" in r.get("_dir", "") or True, "g_100"),]
def rows_for(root, g, filt=lambda r: True):
    return [r for r in load(root) if r["g"] == g and filt(r)]
sets = {
 "101 sp, 50 genes, 100bp (DISCO data)": [r for r in load("/opt/runs/disco") if r["rf"] and "maxgenes" not in r and True],
}
import json, glob
def load_k(root, k):
    out = []
    for f in glob.glob(os.path.join(root, "*-k%d" % k, "*", "result.json")):
        d = json.load(open(f)); d["repid"] = f.split("/")[-2]; out.append(d)
    return out
sets = {
 "101 sp, 50 genes (DISCO data, 5 conditions)": load_k("/opt/runs/disco", 50),
 "101 sp, 1000 genes (DISCO default)": load_k("/opt/runs/disco", 1000),
 "21 sp, 1000 genes, 100bp (QR)": rows_for("/opt/runs/qr", "g_100"),
 "21 sp, 1000 genes, 50bp (QR)": rows_for("/opt/runs/qr", "g_50"),
 "21 sp, 100 genes, 50bp (QR)": rows_for("/opt/runs/qr_k100", "g_50"),
 "21 sp, 1000 true gene trees (QR)": rows_for("/opt/runs/qr", "g_true"),
}
allest = sum((v for k, v in sets.items() if "true" not in k), [])
sets["ALL estimated-gene-tree runs"] = allest
sets["ALL estimated, held-out reps 06-10"] = [r for r in allest if r["repid"] in "06 07 08 09 10".split()]
pairs = [("ASTRID-DISCOR", "ASTRID-DISCO"), ("ASTRAL-DISCOR", "ASTRAL-DISCO"), ("ASTRID-DISCOR-oracle", "ASTRID-DISCO"),
         ("ASTRAL-DISCOR-oracle", "ASTRAL-DISCO"), ("ASTRID-DISCOR-it2", "ASTRID-DISCO"), ("ASTRID-DISCOR-lca", "ASTRID-DISCO"),
         ("ASTRID-DISCOR", "ASTRAL-Pro2"), ("ASTRID-DISCO", "ASTRAL-Pro2")]
print("| experiment | A vs B | n | mean A | mean B | diff | W/T/L | p |\n|---|---|---|---|---|---|---|---|")
for name, R in sets.items():
    for A, B in pairs:
        c = compare(R, A, B)
        if c:
            print(f"| {name} | {A} vs {B} | {c['n']} | {c['mean_a']:.4f} | {c['mean_b']:.4f} | {c['diff']:+.4f} | {c['W']}/{c['T']}/{c['L']} | {c['p']:.3g} |")
