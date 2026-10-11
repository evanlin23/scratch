"""Paired summary of the parsimony pilot (run_blend.py / run_insert.py) and the ML pilot
(run_mlspr_pub.py) on the published conditions.
Usage: python3 pilot_summary.py PARS_DIR ML_DIR OUT_TSV"""
import glob, json, sys, collections
from stats import paired, fmt
rows = []
for f in glob.glob(f"{sys.argv[1]}/*/*/*/blend.json"):
    r = json.load(open(f)); rows.append((r["cond"], r["guide"], r["rep"], "Pars-SPR", r["fn_gtm"], r["fn_blend"], r["constraints_ok"]))
for f in glob.glob(f"{sys.argv[1]}/*/*/*/ins_pars_spr.json"):
    r = json.load(open(f)); rows.append((r["cond"], r["guide"], r["rep"], "Pars-Insert+SPR", r["fn_gtm"], r["fn_spr"], r["ok_spr"]))
for f in glob.glob(f"{sys.argv[2]}/*/*/*/mlspr.json"):
    r = json.load(open(f)); rows.append((r["cond"], r["guide"], r["rep"], "GTM-Blend-ML", r["fn_gtm"], r["fn_mlspr"], r["constraints_ok"]))
with open(sys.argv[3], "w") as o:
    o.write("condition\tguide\treplicate\tmethod\tFN_GTM\tFN_method\tconstraints_ok\n")
    for r in sorted(rows):
        o.write("\t".join(map(str, r)) + "\n")
g = collections.defaultdict(list)
for r in rows:
    g[(r[3], r[0], r[1])].append(r)
for k in sorted(g):
    v = sorted(g[k], key=lambda r: r[2])
    print(f"{k[0]:16s} {k[1]:10s} {k[2]}  A=GTM  " + fmt(paired([r[4] for r in v], [r[5] for r in v])))
for m in sorted({r[3] for r in rows}):
    v = sorted([r for r in rows if r[3] == m])
    print(f"{m:16s} POOLED        A=GTM  " + fmt(paired([r[4] for r in v], [r[5] for r in v])) +
          f"  constraints_ok={all(r[6] for r in v)}")
