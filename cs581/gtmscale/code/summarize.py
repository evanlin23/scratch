"""Collect steps.json (sim, RNASim10K) and result.json (published) into TSVs and paired tests.
Usage: python3 summarize.py RESULTS_DIR
"""
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.stats import wilcoxon

out = sys.argv[1]
rows = []
for p in sorted(glob.glob("/opt/gtms/sim/n*_i*/r*/steps.json") + glob.glob("/opt/gtms/rnasim10k/R*/steps.json")):
    S = json.load(open(p))
    fj = os.path.join(os.path.dirname(p), "full_iq.json")
    if "full_iq" not in S and os.path.exists(fj):
        S["full_iq"] = json.load(open(fj))
    parts = p.split("/")
    cond, rep = ("RNASim10K", parts[-2]) if "/rnasim10k/" in p else (parts[-3], parts[-2])
    guides = sorted({k.split("_m")[0] for k in S if "_m" in k and k.endswith("_subtrees")})
    for g in guides:
        pre = [k[:-len("subtrees")] for k in S if k.startswith(g + "_m") and k.endswith("_subtrees")][0]
        r = dict(cond=cond, rep=rep, guide=g, fn_full_ft=S.get("full_ft", {}).get("fn"),
                 fn_full_iq=S.get("full_iq", {}).get("fn"), t_full_ft=S.get("full_ft", {}).get("wall"),
                 t_full_iq=S.get("full_iq", {}).get("wall"), mem_full_iq=S.get("full_iq", {}).get("mem"),
                 mem_full_ft=S.get("full_ft", {}).get("mem"),
                 fn_guide=S.get("guide_" + g, {}).get("fn"), t_guide=S.get("guide_" + g, {}).get("wall"),
                 fn_subtrees=S[pre + "subtrees"].get("fn"), t_subtrees=S[pre + "subtrees"]["wall"],
                 k=S[pre + "subtrees"]["k"])
        for a in ("gtm", "blendft", "blendfast", "polishft", "cft", "treemerge", "blendml"):
            if pre + a in S and "fn" in S[pre + a]:
                r["fn_" + a] = S[pre + a]["fn"]
                r["t_" + a] = S[pre + a]["wall"]
                r["mem_" + a] = S[pre + a]["mem"]
                if "induced" in S[pre + a]:
                    r["ind_" + a] = "%d/%d" % (S[pre + a]["induced"], S[pre + a]["k"])
        # end-to-end pipeline wall clock (guide + subset trees + GTM [+ blend])
        base = r["t_guide"] or 0
        if g.endswith("+it"):  # add the first round's pipeline cost
            b0 = g[:-3]
            p0 = [k[:-len("subtrees")] for k in S if k.startswith(b0 + "_m") and k.endswith("_subtrees")][0]
            base = (S["guide_" + b0]["wall"] + S[p0 + "subtrees"]["wall"] + S[p0 + "gtm"]["wall"]
                    + S[p0 + "blendft"]["wall"])
        if "t_gtm" in r:
            r["e2e_gtm"] = base + r["t_subtrees"] + r["t_gtm"]
        if "t_blendft" in r:
            r["e2e_blendft"] = r["e2e_gtm"] + r["t_blendft"]
        rows.append(r)

cols = ["cond", "rep", "guide", "k", "fn_guide", "fn_subtrees", "fn_full_ft", "fn_full_iq", "fn_gtm", "fn_blendft",
        "fn_blendfast", "fn_polishft", "fn_cft", "fn_treemerge", "fn_blendml", "ind_blendft", "t_guide", "t_subtrees", "t_gtm", "t_blendft",
        "t_blendfast", "t_polishft", "t_full_ft", "t_full_iq", "e2e_gtm", "e2e_blendft", "mem_blendft", "mem_full_ft", "mem_full_iq"]


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return "%.4f" % v if v < 1.5 else "%.1f" % v
    return str(v)


with open(f"{out}/scale_results.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(fmt(r.get(c)) for c in cols) + "\n")

# published
pub = [json.load(open(p)) for p in sorted(glob.glob("/opt/gtms/pub/*/*/*/result.json"))]
pcols = ["cond", "guide", "rep", "k", "fn_guide", "fn_gtm", "fn_blendft", "fn_polishft", "induced_blendft",
         "wall_blendft", "wall_polishft", "mem_blendft"]
with open(f"{out}/published_results.tsv", "w") as f:
    f.write("\t".join(pcols) + "\n")
    for r in pub:
        f.write("\t".join(fmt(r.get(c)) for c in pcols) + "\n")


def paired(label, a, b, fh):
    """b - a in FN points; negative = b better."""
    x = np.array([(bb - aa) * 100 for aa, bb in zip(a, b) if aa is not None and bb is not None])
    if len(x) == 0:
        return
    w, l_, t = int((x < -1e-9).sum()), int((x > 1e-9).sum()), int((abs(x) <= 1e-9).sum())
    nz = x[abs(x) > 1e-9]
    p = wilcoxon(nz).pvalue if len(nz) >= 1 else float("nan")
    rng = np.random.default_rng(0)
    bs = [rng.choice(x, len(x)).mean() for _ in range(2000)]
    fh.write("%-58s n=%2d mean %+6.2f (95%% CI %+.2f, %+.2f)  better/worse/tie %d/%d/%d  Wilcoxon p=%.3g\n"
             % (label, len(x), x.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5), w, l_, t, p))


with open(f"{out}/paired_tests.txt", "w") as fh:
    fh.write("Differences in FN points (B - A); negative = B more accurate. Two-sided Wilcoxon signed-rank.\n\n")
    groups = {}
    for r in rows:
        groups.setdefault((r["cond"], r["guide"]), []).append(r)
    for (c, g), rs in sorted(groups.items()):
        fh.write("== %s, guide %s (%d reps)\n" % (c, g, len(rs)))
        mean = lambda k: np.mean([r[k] for r in rs if r.get(k) is not None]) * 100 if any(r.get(k) is not None for r in rs) else float("nan")
        fh.write("   mean FN%%: guide %.2f subsets %.2f fullFT %.2f fullIQ %.2f GTM %.2f BlendFT %.2f polishFT %.2f cFT %.2f TreeMerge %.2f\n"
                 % tuple(mean(k) for k in ("fn_guide", "fn_subtrees", "fn_full_ft", "fn_full_iq", "fn_gtm", "fn_blendft",
                                           "fn_polishft", "fn_cft", "fn_treemerge")))
        tm = lambda k: np.mean([r[k] for r in rs if r.get(k) is not None]) if any(r.get(k) is not None for r in rs) else float("nan")
        fh.write("   mean wall s: e2e GTM %.0f, e2e BlendFT %.0f, full FT %.0f, full IQ(-fast) %.0f\n"
                 % (tm("e2e_gtm"), tm("e2e_blendft"), tm("t_full_ft"), tm("t_full_iq")))
        fh.write("   mean FN%%: BlendFast %.2f; mean wall s BlendFT %.0f BlendFast %.0f\n" % (mean("fn_blendfast"), tm("t_blendft"), tm("t_blendfast")))
        for A, B in (("gtm", "blendft"), ("gtm", "blendfast"), ("gtm", "polishft"), ("polishft", "blendft"), ("gtm", "treemerge"),
                     ("full_ft", "blendft"), ("full_iq", "blendft"), ("full_ft", "gtm")):
            paired("   %s -> %s" % (A, B), [r.get("fn_" + A) for r in rs], [r.get("fn_" + B) for r in rs], fh)
    # pooled over simulated ftfast conditions
    fh.write("\n== published conditions (published guide + IQ-TREE subset trees + published GTM tree)\n")
    pg = {}
    for r in pub:
        pg.setdefault((r["cond"], r["guide"]), []).append(r)
    for (c, g), rs in sorted(pg.items()):
        fh.write("== %s, guide %s (%d reps): GTM %.2f BlendFT %.2f polishFT %.2f\n" % (
            c, g, len(rs), *(np.mean([r[k] for r in rs]) * 100 for k in ("fn_gtm", "fn_blendft", "fn_polishft"))))
        for A, B in (("gtm", "blendft"), ("gtm", "polishft"), ("polishft", "blendft")):
            paired("   %s -> %s" % (A, B), [r["fn_" + A] for r in rs], [r["fn_" + B] for r in rs], fh)
    fh.write("\n== POOLED (each replicate x guide x condition is one pair)\n")
    simrows = [r for r in rows if r["cond"] != "RNASim10K"]
    first = [r for r in simrows if not r["guide"].endswith("+it")]
    for lab, rs in (("simulated, first round (ftfast + kmer guides)", first),
                    ("simulated, iteration round", [r for r in simrows if r["guide"].endswith("+it")]),
                    ("published (1000M1-HF, Cox1-HET; FT + IQ guides)", pub)):
        for A, B in (("gtm", "blendft"), ("gtm", "blendfast"), ("gtm", "polishft"), ("polishft", "blendft")):
            paired("   %s: %s -> %s" % (lab, A, B), [r.get("fn_" + A) for r in rs], [r.get("fn_" + B) for r in rs], fh)
    allp = first + pub
    paired("   ALL first-round + published: gtm -> blendft", [r.get("fn_gtm") for r in allp],
           [r.get("fn_blendft") for r in allp], fh)
print(open(f"{out}/paired_tests.txt").read())
