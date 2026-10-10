"""Tables for REPORT.md.  Usage: python analyze.py > ../results/tables.md"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
MAIN = ["astrid-pro", "astrid-disco", "asteroid", "astral-pro3", "wqfm-gdl"]
KEY = ["data", "cond", "rep", "ngen", "level"]


def load(*files):
    rows = []
    for f in files:
        p = os.path.join(R, f)
        if not os.path.exists(p):
            continue
        for l in open(p):
            r = json.loads(l)
            r["src"] = f
            rows.append(r)
    d = pd.DataFrame(rows)
    for k in KEY:
        if k not in d:
            d[k] = None
    d["ngen"] = d["ngen"].astype(str)
    return d


def ok(d):
    return d[d["error"].isna()] if "error" in d else d


def wide(d, idx=KEY):
    d = ok(d).drop_duplicates(idx + ["method"], keep="last")
    return d.pivot_table(index=idx, columns="method", values="FNrate", aggfunc="first")


def pair(w, a, b):
    x = w[[a, b]].dropna() if a in w and b in w else pd.DataFrame(columns=[a, b])
    n = len(x)
    if n == 0:
        return None
    dd = x[a] - x[b]
    W, T, L = int((dd < 0).sum()), int((dd == 0).sum()), int((dd > 0).sum())
    try:
        p = wilcoxon(dd[dd != 0]).pvalue if (dd != 0).sum() > 0 else 1.0
    except ValueError:
        p = 1.0
    return dict(n=n, diff=dd.mean(), W=W, T=T, L=L, p=p)


def holm(ps):
    ps = np.array(ps, dtype=float)
    o = np.argsort(ps)
    m = len(ps)
    adj = np.empty(m)
    run = 0
    for r, i in enumerate(o):
        run = max(run, (m - r) * ps[i])
        adj[i] = min(1, run)
    return adj


def fmt(x, nd=3):
    return "–" if x is None or (isinstance(x, float) and np.isnan(x)) else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def table(head, rows):
    out = "| " + " | ".join(head) + " |\n|" + "---|" * len(head) + "\n"
    for r in rows:
        out += "| " + " | ".join(fmt(x) if not isinstance(x, str) else x for x in r) + " |\n"
    return out


def means(w, groups, cols):
    rows = []
    for g, sub in w.groupby(level=groups):
        g = g if isinstance(g, tuple) else (g,)
        best = min((sub[c].mean() for c in cols if c in sub and sub[c].notna().any()), default=np.nan)
        r = list(g) + [len(sub)]
        for c in cols:
            if c in sub and sub[c].notna().any():
                m, k = sub[c].mean(), sub[c].notna().sum()
                s = f"{m:.3f}" + ("" if k == len(sub) else f" ({k})")
                r.append(f"**{s}**" if abs(m - best) < 1e-9 else s)
            else:
                r.append("–")
        rows.append(r)
    return table(list(groups) + ["n"] + cols, rows)


def comparisons(w, ref, others, title):
    res = [(o, pair(w, ref, o)) for o in others]
    res = [(o, r) for o, r in res if r]
    adj = holm([r["p"] for _, r in res]) if res else []
    rows = [[o, r["n"], f"{r['diff']:+.4f}", f"{r['W']}/{r['T']}/{r['L']}", f"{r['p']:.2g}", f"{a:.2g}"]
            for (o, r), a in zip(res, adj)]
    return f"**{title}**\n\n" + table([f"{ref} vs", "n", "mean diff", "W/T/L", "p", "Holm p"], rows)


def oracle(w, cols):
    cols = [c for c in cols if c in w]
    x = w[cols].dropna()
    if len(x) == 0:
        return None
    best_single = x.mean().min()
    return dict(n=len(x), best=x.mean().idxmin(), best_single=best_single, oracle=x.min(axis=1).mean(),
                worst=x.mean().max())


def main():
    out = sys.stdout
    # ---------------- Q1 decomposition, FastMulRFS
    d = load("fmrfs.jsonl")
    if len(d):
        w = wide(d)
        out.write("## Q1. Error decomposition: FastMulRFS data (100 taxa, 100 genes, 6 conditions x 10 reps)\n\n")
        lv = w.reset_index()
        rows = []
        for m in ["astrid-pro", "astrid-disco", "disco-astral", "astrid-multi", "asteroid", "astral-pro3", "wqfm-gdl"]:
            r = [m]
            tt = {"astrid-pro": "astrid-pro-tt", "astrid-disco": "astrid-disco-tt", "disco-astral": "disco-astral-tt"}.get(m)
            g = lv.groupby("level")
            L0 = g[tt].mean().get("true", np.nan) if tt and tt in lv else np.nan
            L1 = g[m].mean().get("true", np.nan) if m in lv else np.nan
            L2a = g[m].mean().get("est100", np.nan) if m in lv else np.nan
            L2b = g[m].mean().get("est25", np.nan) if m in lv else np.nan
            r += [L0, L1, L2a, L2b, L1 - L0 if tt else np.nan, L2a - L1, L2b - L1]
            rows.append(r)
        out.write(table(["method", "L0 true tree+tags", "L1 true tree", "L2 est 100bp", "L2 est 25bp",
                         "tagging (L1-L0)", "GTEE 100bp (L2-L1)", "GTEE 25bp"], rows) + "\n")
        out.write("Mean FN rate by condition and level:\n\n")
        out.write(means(w, ["level", "cond"], ["astrid-pro-tt", "astrid-pro", "astrid-disco-tt", "astrid-disco",
                                               "disco-astral-tt", "disco-astral", "astrid-multi", "asteroid",
                                               "astral-pro3", "wqfm-gdl"]) + "\n")
        for lev in ["true", "est100", "est25"]:
            sub = w.xs(lev, level="level")
            out.write(comparisons(sub, "astrid-pro", MAIN[1:] + ["astrid-multi", "disco-astral"],
                                  f"FastMulRFS level {lev}: ASTRID-Pro minus method") + "\n")
        # level effects per method
        rows = []
        for m in MAIN + ["astrid-multi", "disco-astral"]:
            for a, b in [("est100", "true"), ("est25", "est100")]:
                x = lv.pivot_table(index=["cond", "rep"], columns="level", values=m)
                if a in x and b in x:
                    p = pair(x, a, b)
                    rows.append([m, f"{a} - {b}", p["n"], f"{p['diff']:+.4f}", f"{p['W']}/{p['T']}/{p['L']}", f"{p['p']:.2g}"])
            tt = m + "-tt"
            if tt in lv:
                x = lv[lv.level == "true"].set_index(["cond", "rep"])
                p = pair(x, m, tt)
                rows.append([m, "true tags: est - true", p["n"], f"{p['diff']:+.4f}", f"{p['W']}/{p['T']}/{p['L']}", f"{p['p']:.2g}"])
        out.write("**Level effects (W = first level has lower error)**\n\n" +
                  table(["method", "contrast", "n", "mean diff", "W/T/L", "p"], rows) + "\n")
        rows = []
        for lev in ["true", "est100", "est25"]:
            sub = w.xs(lev, level="level")
            o = oracle(sub, MAIN)
            o2 = oracle(sub, MAIN + ["astrid-multi", "disco-astral"])
            if o:
                rows.append([lev, o["n"], o["best"], o["best_single"], o["oracle"], o["best_single"] - o["oracle"],
                             o["worst"] - o["best_single"], o2["oracle"] if o2 else np.nan])
        out.write("**Headroom: best single method vs per-replicate oracle over the 5 main methods**\n\n" +
                  table(["level", "n", "best method", "best mean FN", "oracle mean FN", "oracle gain",
                         "worst - best", "oracle over 7"], rows) + "\n")
        # hybrid
        rows = []
        for lev in ["est100", "est25"]:
            sub = w.xs(lev, level="level")
            tm = ok(d[d.level == lev]).groupby("method")["sec"].mean()
            for m in ["apro3-fast", "apro3-guide-fast", "apro3-guide", "astrid-pro"]:
                p = pair(sub, m, "astral-pro3")
                if p:
                    rows.append([lev, m, p["n"], sub[m].mean(), sub["astral-pro3"].mean(), f"{p['diff']:+.4f}",
                                 f"{p['W']}/{p['T']}/{p['L']}", f"{p['p']:.2g}", tm.get(m, np.nan), tm.get("astral-pro3")])
        out.write("**Q4 hybrid on FastMulRFS (diff = hybrid - ASTRAL-Pro3; seconds = mean wall time, 1 thread)**\n\n" +
                  table(["level", "variant", "n", "FN variant", "FN ASTRAL-Pro3", "diff", "W/T/L", "p", "sec variant",
                         "sec ASTRAL-Pro3"], rows) + "\n")
    # ---------------- Q1 decomposition, DISCO
    d = load("disco_q1.jsonl")
    if len(d):
        prior = load("prior/disco_runs.jsonl")
        prior = prior[prior.method.isin(["astral-pro3", "astrid-pro", "astrid-disco", "astrid-multi", "asteroid"])].copy()
        prior["level"] = "est100"
        prior["ngen"] = "1000"
        d = pd.concat([prior[prior.cond.isin(d.cond.unique()) & prior.rep.isin(d.rep.unique())], d])
        w = wide(d)
        out.write("## Q1. Error decomposition: DISCO data (100 taxa, 1000 genes)\n\n")
        out.write("Estimated-tree (est100) results not rerun here are taken from the previous pilot (same software "
                  "versions, same inputs; ASTRAL-Pro3 deterministic seed).\n\n")
        out.write(means(w, ["cond", "level"], ["astrid-pro-tt", "astrid-pro", "astrid-disco-tt", "astrid-disco",
                                               "astrid-multi", "asteroid", "astral-pro3", "wqfm-gdl"]) + "\n")
        for lev in ["true", "est100"]:
            sub = w.xs(lev, level="level")
            out.write(comparisons(sub, "astrid-pro", MAIN[1:] + ["astrid-multi"], f"DISCO level {lev}") + "\n")
            o = oracle(sub, MAIN)
            if o:
                out.write(f"Oracle over 5 methods ({lev}, n={o['n']}): best single {o['best']} {o['best_single']:.4f}, "
                          f"oracle {o['oracle']:.4f}.\n\n")
    d = load("disco_true.jsonl", "disco_q1.jsonl")
    if len(d):
        d = d[d.level == "true"]
        prior = load("prior/disco_runs.jsonl")
        prior["level"] = "est100"
        prior["ngen"] = "1000"
        d["ngen"] = "1000"
        d = pd.concat([prior, d])
        w = wide(d)
        cols = ["astrid-pro", "astrid-disco", "astrid-multi", "asteroid"]
        rows = []
        lv = w.reset_index()
        for c, g in lv.groupby("cond"):
            x = g.pivot_table(index="rep", columns="level", values=cols)
            r = [c]
            for m in cols:
                if (m, "true") in x and (m, "est100") in x:
                    y = x[m][["true", "est100"]].dropna()
                    r += [f"{y['true'].mean():.3f} → {y['est100'].mean():.3f} ({len(y)})"]
                else:
                    r += ["–"]
            rows.append(r)
        out.write("## Q1. DISCO: true-tree (L1) → estimated-tree (L2, 100 bp) FN rate per condition, fast methods, "
                  "reps 01-07 (n in parentheses)\n\n" + table(["cond"] + cols, rows) + "\n")
        sub = w.xs("true", level="level")
        out.write(comparisons(sub, "astrid-pro", ["astrid-disco", "asteroid", "astrid-multi"], "DISCO true trees, all conditions") + "\n")
    # ---------------- tag accuracy
    for f in ["tagacc_fmrfs.jsonl", "tagacc_disco.jsonl"]:
        p = os.path.join(R, f)
        if not os.path.exists(p):
            continue
        rows = []
        recs = [json.loads(l) for l in open(p)]
        grp = {}
        for r in recs:
            c = r["file"].split("/")[-3]
            grp.setdefault(c, []).append(r)
        for c, rs in sorted(grp.items()):
            row = [c, len(rs)]
            for m in ["astral-pro3", "disco"]:
                first = m == "astral-pro3"
                tp = sum(r[m]["true_para_pairs"] for r in rs)
                fo = sum(r[m]["false_orth_rate"] * r[m]["true_para_pairs"] for r in rs)
                to = sum(r[m]["true_orth_pairs"] for r in rs)
                fp = sum(r[m]["false_para_rate"] * r[m]["true_orth_pairs"] for r in rs)
                row += ([tp / (tp + to)] if first else []) + [fo / max(tp, 1), fp / max(to, 1)]
            rows.append(row)
        out.write(f"## Tag accuracy on true gene trees ({f})\n\n" + table(
            ["cond", "reps", "paralog pair frac", "A-Pro3 false-orth", "A-Pro3 false-para", "DISCO/MinDup false-orth",
             "DISCO/MinDup false-para"], rows) + "\n")
    # ---------------- Q2 regimes
    d = load("prior/disco_runs.jsonl", "disco_q2.jsonl")
    if len(d):
        d["level"] = d["level"].fillna("est100")
        d["ngen"] = "1000"
        w = wide(d)
        cols = ["astrid-pro", "astrid-multi", "astrid-disco", "asteroid", "astral-pro3", "wqfm-gdl", "disco-astral"]
        out.write("## Q2. DISCO regimes (estimated gene trees, 100 bp, 1000 genes)\n\n")
        out.write(means(w, ["cond"], cols) + "\n")
        rows = []
        for c, sub in w.groupby(level="cond"):
            for o in ["astrid-disco", "asteroid", "astral-pro3", "wqfm-gdl", "astrid-multi"]:
                p = pair(sub, "astrid-pro", o)
                if p:
                    rows.append([c, o, p["n"], f"{p['diff']:+.4f}", f"{p['W']}/{p['T']}/{p['L']}", f"{p['p']:.2g}"])
        out.write("**ASTRID-Pro minus method, per condition**\n\n" + table(["cond", "vs", "n", "mean diff", "W/T/L", "p"], rows) + "\n")
    d = load("genes_q2.jsonl")
    if len(d):
        w = wide(d)
        out.write("## Q2. Few vs many genes (DISCO gtrees_10000_l1, 100 taxa, 100 bp)\n\n")
        out.write(means(w, ["level", "ngen"], ["astrid-pro", "astrid-multi", "astrid-disco", "asteroid", "astral-pro3", "wqfm-gdl"]) + "\n")
        rows = []
        for g, sub in w.xs("est100", level="level").groupby(level="ngen"):
            for o in ["astrid-disco", "asteroid", "astral-pro3", "wqfm-gdl"]:
                p = pair(sub, "astrid-pro", o)
                if p:
                    rows.append([g, o, p["n"], f"{p['diff']:+.4f}", f"{p['W']}/{p['T']}/{p['L']}", f"{p['p']:.2g}"])
        out.write(table(["ngen", "ASTRID-Pro vs", "n", "mean diff", "W/T/L", "p"], rows) + "\n")
    d = load("sp1000_q2.jsonl")
    if len(d):
        w = wide(d)
        out.write("## Q2. 1000 species (DISCO species_1000)\n\n")
        out.write(means(w, ["level", "ngen"], ["astrid-pro-tt", "astrid-pro", "astrid-multi", "astrid-disco", "asteroid",
                                              "astral-pro3", "wqfm-gdl"]) + "\n")
        for lev in sorted(set(w.index.get_level_values("level"))):
            sub = w.xs(lev, level="level")
            out.write(comparisons(sub, "astrid-pro", ["astrid-disco", "asteroid", "astrid-multi", "astral-pro3"],
                                  f"species_1000, level {lev}") + "\n")
    d = load("hybrid_disco.jsonl")
    if len(d):
        prior = load("prior/disco_runs.jsonl")
        prior = prior[prior.method.isin(["astral-pro3", "astrid-pro"])].copy()
        prior["level"] = "est100"
        prior["ngen"] = "1000"
        d = pd.concat([d, prior[prior.cond.isin(d.cond.unique()) & prior.rep.isin(d.rep.unique())]])
        w = wide(d)
        tm = ok(d).groupby("method")["sec"].mean()
        rows = []
        for m in ["apro3-fast", "apro3-guide-fast", "astrid-pro"]:
            p = pair(w, m, "astral-pro3")
            if p:
                rows.append([m, p["n"], w[m].mean(), f"{p['diff']:+.4f}", f"{p['W']}/{p['T']}/{p['L']}", f"{p['p']:.2g}",
                             tm.get(m, np.nan)])
        out.write("## Q4. Hybrid on DISCO (1000 genes, 100 bp; ASTRAL-Pro3 from previous pilot, mean "
                  f"{tm.get('astral-pro3', np.nan):.0f} s)\n\n" +
                  table(["variant", "n", "FN", "diff vs ASTRAL-Pro3", "W/T/L", "p", "sec"], rows) + "\n")
    d = load("scaling.jsonl")
    if len(d):
        d = d.copy()
        d["cfg"] = d["method"] + "@" + d["threads"].astype(str) + "t"
        d["cond"] = d["cond"].where(d["cond"] != "genes", "genes")
        dd = d.drop_duplicates(["cond", "ngen", "cfg"], keep="last")
        out.write("## Q3. Scaling (idle machine, one job at a time)\n\n")
        for v, lab in [("sec", "wall time (s)"), ("rssMB", "peak RSS (MB)")]:
            t = dd.pivot_table(index=["cond", "ngen"], columns="cfg", values=v, aggfunc="first")
            if "error" in dd:
                e = dd[dd["error"].notna()]
                for _, r in e.iterrows():
                    t.loc[(r["cond"], r["ngen"]), r["cfg"]] = np.nan
            out.write(f"**{lab}**\n\n" + t.round(1).fillna("–").to_markdown() + "\n\n")
        if "FNrate" in dd:
            t = dd.pivot_table(index=["cond", "ngen"], columns="cfg", values="FNrate", aggfunc="first")
            out.write("**FN rate**\n\n" + t.round(3).fillna("–").to_markdown() + "\n\n")
    p = os.path.join(R, "qscore.jsonl")
    if os.path.exists(p):
        q = pd.DataFrame([json.loads(l) for l in open(p)]).drop_duplicates(["cond", "rep", "level"], keep="last")
        rows = []
        for lev, g in q.groupby("level"):
            a, t, b = g["astral-pro3_score"], g["true_score"], g["astrid-pro_score"]
            rows.append([lev, len(g), int((t > a).sum()), int((t == a).sum()), int((t < a).sum()),
                         float(((a - t) / a).mean() * 100), int((b < a).sum()),
                         int(((g["astrid-pro_FN"] < g["astral-pro3_FN"]) & (b < a)).sum()),
                         float(g["astral-pro3_FN"].mean()), float(g["astrid-pro_FN"].mean())])
        out.write("## Q4b. Search or objective? ASTRAL-Pro3 quartet score of the true tree vs ASTRAL-Pro3's tree "
                  "(FastMulRFS, 100 genes, reps 01-05)\n\n" + table(
                      ["level", "n", "true scores higher (search failure)", "equal", "true scores lower (objective failure)",
                       "score gap ASTRAL - true (%)", "ASTRID-Pro tree scores lower than ASTRAL's",
                       "... and has fewer FN", "mean FN ASTRAL-Pro3", "mean FN ASTRID-Pro"], rows) + "\n")
    # failures
    rows = []
    for f in ["fmrfs.jsonl", "disco_q1.jsonl", "disco_q2.jsonl", "genes_q2.jsonl", "sp1000_q2.jsonl", "hybrid_disco.jsonl",
              "scaling.jsonl"]:
        d = load(f)
        if len(d) and "error" in d:
            e = d[d["error"].notna()]
            for (m, c, er), g in e.groupby(["method", "cond", "error"]):
                rows.append([f, m, c, er, len(g)])
    if rows:
        out.write("## Failed runs\n\n" + table(["file", "method", "cond", "error", "n"], rows) + "\n")


if __name__ == "__main__":
    main()
