"""Analyse md_exp.jsonl. Paired by dataset (condition x replicate x #genes x gene-tree type).

Split: DEV = odd replicates (used to pick the variant), TEST = even replicates (held out).
For each comparison method vs baseline: mean FN difference (method - baseline; negative = better),
W/T/L counted per dataset with tie band |diff| < 1e-9 (FN is a multiple of 1/23, so ties are
exact equality), two-sided Wilcoxon signed-rank p (zero differences dropped, 'wilcox').
"""
import os, json, sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

RES = os.path.join(os.path.dirname(__file__), "..", "results")


def compare(d, a, b, keys):
    x = d[d.method == a].set_index(keys).FN
    y = d[d.method == b].set_index(keys).FN
    j = x.index.intersection(y.index)
    diff = (x.loc[j] - y.loc[j]).values
    if len(diff) == 0:
        return None
    w = int((diff < -1e-9).sum()); l = int((diff > 1e-9).sum()); t = len(diff) - w - l
    p = wilcoxon(diff[np.abs(diff) > 1e-9]).pvalue if (w + l) > 0 else 1.0
    return dict(method=a, vs=b, n=len(diff), mean_FN_a=x.loc[j].mean(), mean_FN_b=y.loc[j].mean(),
                mean_diff=diff.mean(), W=w, T=t, L=l, p=p)


def main():
    d = pd.DataFrame([json.loads(l) for l in open("/opt/runs/md_exp.jsonl")])
    d["split"] = np.where(d.rep % 2 == 1, "dev", "test")
    keys = ["height", "rate", "pattern", "ngenes", "rep", "genetrees"]
    rows = []
    for split in ["dev", "test"]:
        for pat in sorted(d.pattern.unique()) + ["ALL"]:
            for gt in ["true", "raxml"]:
                s = d[(d.split == split) & (d.genetrees == gt)]
                if pat != "ALL":
                    s = s[s.pattern == pat]
                for a, b in [("astrid-ht", "astrid"), ("astrid-plugin", "astrid"), ("astrid-cmp", "astrid"),
                             ("asteroid", "astrid"), ("astral", "astrid"), ("astrid-plugin", "asteroid"),
                             ("astrid-plugin", "astral"), ("astrid-cmp", "astral"), ("astrid", "astrid-pub")]:
                    r = compare(s, a, b, keys)
                    if r:
                        r.update(split=split, pattern=pat, genetrees=gt)
                        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/md_wilcoxon.csv", index=False)
    # mean FN table
    tab = d.pivot_table(index=["split", "pattern", "genetrees", "ngenes"], columns="method", values="FN", aggfunc="mean")
    tab.round(4).to_csv(f"{RES}/md_meanFN.csv")
    rt = d.groupby(["method", "ngenes"]).sec.mean().unstack()
    rt.round(3).to_csv(f"{RES}/md_runtime.csv")
    pd.set_option("display.width", 250)
    print(d.groupby("split").rep.nunique(), len(d))
    print(tab.round(3).to_string())
    print(out[out.pattern == "ALL"].round(4).to_string())
    print(rt.round(2))


if __name__ == "__main__":
    main()
