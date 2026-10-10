"""Score every run against its BAliBASE reference and print paired comparisons.

Usage: python3 summarize.py RUNROOT OUT.csv [baseline]
"""
import glob, os, sys, csv
import numpy as np
from scipy.stats import wilcoxon
from bbscore import score

root, out = sys.argv[1], sys.argv[2]
base = sys.argv[3] if len(sys.argv) > 3 else "linsi"
# oracle (experimental-structure 3Di) rows only where every sequence mapped to a PDB chain (min coverage >= 0.9)
cov = {}
covf = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "true3di_coverage.tsv")
if os.path.exists(covf):
    for line in open(covf).read().splitlines()[1:]:
        f = line.split("\t"); cov[f[0]] = float(f[2])
rows = []
for d in sorted(glob.glob(os.path.join(root, "BB*"))):
    sid = os.path.basename(d)
    grp = "RV" + sid[-5:-3] + ("-BBS" if sid.startswith("BBS") else "-BB")
    xml = f"/opt/bb3/bb3_release/RV{sid[-5:-3]}/{sid}.xml"
    for f in sorted(glob.glob(os.path.join(d, "*.fa"))):
        m = os.path.basename(f)[:-3]
        if m.endswith(".3di") or (m.endswith("_true") and cov.get(sid, 0) < 0.9):
            continue
        try:
            sp, tc, fn, fp = score(xml, f)
        except BaseException as e:
            print("FAIL", sid, m, e, file=sys.stderr); continue
        t = open(os.path.join(d, m + ".time")).read().split()[-3:] if os.path.exists(os.path.join(d, m + ".time")) else ["nan"] * 3
        rows.append(dict(set=sid, group=grp, method=m, SP=sp, TC=tc, SPFN=fn, SPFP=fp,
                         wall=float(t[0]), cpu=float(t[1]) + float(t[2])))
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

import collections
tab = collections.defaultdict(dict)
for r in rows:
    tab[(r["group"], r["set"])][r["method"]] = r
methods = sorted({r["method"] for r in rows})
groups = sorted({r["group"] for r in rows}) + sorted({r["group"][:4] for r in rows}) + ["ALL"]
for g in groups:
    keys = [k for k in tab if g == "ALL" or k[0] == g or k[0].startswith(g + "-")]
    print(f"\n== {g} ==")
    print(f"{'method':10s} {'n':>3s} {'SP':>6s} {'TC':>6s} | vs {base}: dSP  W/T/L  p(SP)   dTC  W/T/L  p(TC)")
    for m in methods:
        ks = [k for k in keys if m in tab[k] and base in tab[k]]
        if not ks: continue
        sp = np.array([tab[k][m]["SP"] for k in ks]); tc = np.array([tab[k][m]["TC"] for k in ks])
        bsp = np.array([tab[k][base]["SP"] for k in ks]); btc = np.array([tab[k][base]["TC"] for k in ks])
        def wtl(a, b):
            d = a - b; return f"{(d>1e-9).sum()}/{(abs(d)<=1e-9).sum()}/{(d<-1e-9).sum()}"
        def p(a, b):
            d = a - b
            return wilcoxon(d).pvalue if (abs(d) > 1e-9).sum() >= 5 else float("nan")
        print(f"{m:10s} {len(ks):3d} {sp.mean():6.3f} {tc.mean():6.3f} | {sp.mean()-bsp.mean():+.3f} {wtl(sp,bsp):>8s} {p(sp,bsp):.1e} "
              f"{tc.mean()-btc.mean():+.3f} {wtl(tc,btc):>8s} {p(tc,btc):.1e}")
