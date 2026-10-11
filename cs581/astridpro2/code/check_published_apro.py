"""Sanity check of our ASTRAL-Pro3 runs: score the ASTRAL-Pro (v2) species trees published by the
ASTRAL-Pro authors (github.com/chaoszhang/A-pro_data, S100 = FastMulRFS data) on the same replicates,
sequence lengths and gene counts, and compare with our astral-pro3 and astrid-pro FN rates.
Usage: python check_published_apro.py results/fmrfs_runs.jsonl"""
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "gdl", "code"))
from phylo import parse_newick, rf_error
import numpy as np
from scipy.stats import wilcoxon
P = "/opt/data/apro_pub"
ours = {}
for l in open(sys.argv[1]):
    r = json.loads(l)
    if "error" in r:
        continue
    ours[(r["cond"], r["rep"], r["sqln"], r["ngen"], r["method"])] = r["FNrate"]
rows = []
for (cond, rep, sq, ng, m), v in ours.items():
    if m != "astral-pro3":
        continue
    dl = {"1": "0.0000000001", "2": "0.0000000002", "5": "0.0000000005"}[cond[2]]
    ps = "10000000" if cond.endswith("1e7") else "50000000"
    d = f"{P}/ntaxa-100.dlrate-{dl}.psize-{ps}/{rep}"
    f = f"{d}/apro-v2-raxml-sqln-{sq}-ngen-{ng}.tree"
    if not os.path.exists(f):
        continue
    true = parse_newick(open(f"/opt/data/fmrfs/ntaxa-100.dlrate-{dl}.psize-{ps}/{rep}/s_tree.trees").read().strip())
    est = parse_newick(open(f).read().strip().split("\n")[0])
    fn, fp, i1, i2 = rf_error(est, true)
    rows.append((fn / i1, v, ours.get((cond, rep, sq, ng, "astrid-pro"))))
a = np.array(rows, dtype=float)
print(f"n={len(a)}  published ASTRAL-Pro v2: {a[:,0].mean():.4f}  our ASTRAL-Pro3: {a[:,1].mean():.4f}  "
      f"our ASTRID-Pro: {np.nanmean(a[:,2]):.4f}")
d = a[:, 2] - a[:, 0]
print(f"ASTRID-Pro - published ASTRAL-Pro: mean {d.mean():+.4f}, W/T/L {(d<0).sum()}/{(d==0).sum()}/{(d>0).sum()}, "
      f"p={wilcoxon(d[d!=0]).pvalue:.2g}")
d2 = a[:, 1] - a[:, 0]
print(f"our ASTRAL-Pro3 - published ASTRAL-Pro: mean {d2.mean():+.4f}, p={wilcoxon(d2[d2!=0]).pvalue:.2g}")
