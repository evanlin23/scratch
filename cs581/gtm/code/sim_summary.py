"""Aggregate simulation results (mlspr.json per replicate) with paired statistics.
Usage: python3 sim_summary.py SIMROOT OUT_TSV"""
import glob, json, os, sys
from stats import paired, fmt
rows = []
for f in sorted(glob.glob(f"{sys.argv[1]}/*/*/sub_*/mlspr.json")):
    r = json.load(open(f))
    sd = os.path.dirname(f)
    if os.path.exists(f"{sd}/treemerge.json"):
        r["fn_treemerge"] = json.load(open(f"{sd}/treemerge.json"))["fn_treemerge"]
    rd = os.path.dirname(sd)
    if os.path.exists(f"{rd}/iqtree_full.json"):
        r["fn_iqtree_full"] = json.load(open(f"{rd}/iqtree_full.json"))["fn_iqtree_full"]
    cond = f.split("/")[-4] + "/" + r["subtrees"]
    rows.append((cond, r))
with open(sys.argv[2], "w") as o:
    o.write("condition\trep\tk\tFN_guide\tFN_subsets\tFN_GTM\tFN_best_unblended\tFN_oracle_blend\tFN_GTM_Blend_ML\t"
            "FN_TreeMerge\tFN_IQTREE_full\tmoves\tlogL_GTM\tlogL_GTM_Blend_ML\tlogL_true\tconstraints_ok\tseconds\n")
    for c, r in rows:
        o.write("\t".join(str(x) for x in [c, r["rep"], r.get("k"), r.get("fn_guide"), r.get("fn_subsets"), r["fn_gtm"],
                                          r.get("fn_best_unblended"), r.get("fn_oracle_ins"), r["fn_mlspr"], r.get("fn_treemerge"), r.get("fn_iqtree_full"), r["moves"],
                                          r["logL_gtm"], r["logL_mlspr"], r["logL_true"], r["constraints_ok"],
                                          round(r["sec"])]) + "\n")
for c in sorted({c for c, _ in rows}):
    rs = sorted([r for cc, r in rows if cc == c], key=lambda r: r["rep"])
    print(f"== {c}  (n={len(rs)}, all constraints ok: {all(r['constraints_ok'] for r in rs)}, "
          f"mean guide FN {100*sum(r['fn_guide'] for r in rs)/len(rs):.1f}%, subset FN {100*sum(r['fn_subsets'] for r in rs)/len(rs):.1f}%, "
          f"best-unblended {100*sum(r['fn_best_unblended'] for r in rs)/len(rs):.1f}%, oracle-blend {100*sum(r['fn_oracle_ins'] for r in rs)/len(rs):.1f}%, "
          f"mean sec {sum(r['sec'] for r in rs)/len(rs):.0f})")
    print("  FN   A=GTM B=GTM-Blend-ML  " + fmt(paired([r["fn_gtm"] for r in rs], [r["fn_mlspr"] for r in rs])))
    tm = [r for r in rs if "fn_treemerge" in r]
    if tm:
        print("  FN   A=GTM B=TreeMerge     " + fmt(paired([r["fn_gtm"] for r in tm], [r["fn_treemerge"] for r in tm])))
        print("  FN   A=TreeMerge B=GTM-Blend-ML " + fmt(paired([r["fn_treemerge"] for r in tm], [r["fn_mlspr"] for r in tm])))
    iq = [r for r in rs if "fn_iqtree_full" in r]
    if iq:
        print("  FN   A=IQ-TREE(full data) B=GTM-Blend-ML " + fmt(paired([r["fn_iqtree_full"] for r in iq], [r["fn_mlspr"] for r in iq])))
    print("  logL improvement mean %.1f; logL(true) - logL(Blend-ML) mean %.1f" % (
        sum(r["logL_mlspr"] - r["logL_gtm"] for r in rs) / len(rs), sum(r["logL_true"] - r["logL_mlspr"] for r in rs) / len(rs)))
