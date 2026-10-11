"""Tables from bench.py results.
Usage: python analyze.py RESULTS.jsonl REF_METHOD [--threads 1] [--md OUT.md] [--title T] [--holm m1,m2,...]
Prints (and optionally writes) mean FN rate per method (overall, per data set, per condition x sqln), runtime, and
paired comparisons REF - other (negative = REF better): n, mean diff, 95% bootstrap CI, W/T/L, two-sided
Wilcoxon p (zero differences dropped), Holm-adjusted p over the methods listed in --holm."""
import argparse
import collections
import json
import random

from scipy.stats import wilcoxon

KEY = ("data", "cond", "rep", "sqln")


def load(path, threads):
    R = collections.defaultdict(dict)
    T = collections.defaultdict(dict)
    err = collections.Counter()
    for l in open(path):
        r = json.loads(l)
        if r.get("threads", 1) != threads:
            continue
        k = tuple(r[x] for x in KEY)
        if "error" in r:
            err[(r["method"], r["data"], r["cond"], r["error"])] += 1
            continue
        R[r["method"]][k] = r["FNrate"]
        T[r["method"]][k] = (r["wall"], r["cpu"])
    return R, T, err


def boot(d, B=4000, seed=1):
    rnd = random.Random(seed)
    n = len(d)
    ms = sorted(sum(d[rnd.randrange(n)] for _ in range(n)) / n for _ in range(B))
    return ms[int(0.025 * B)], ms[int(0.975 * B) - 1]


def compare(R, ref, other):
    ks = sorted(set(R[ref]) & set(R[other]))
    d = [R[ref][k] - R[other][k] for k in ks]
    if not d:
        return None
    w = sum(x < -1e-12 for x in d)
    l_ = sum(x > 1e-12 for x in d)
    t = len(d) - w - l_
    p = 1.0
    if w + l_ > 0:
        p = wilcoxon(d, zero_method="wilcox").pvalue
    lo, hi = boot(d)
    return dict(n=len(d), mean=sum(d) / len(d), lo=lo, hi=hi, W=w, T=t, L=l_, p=p)


def holm(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj = [0.0] * len(ps)
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - r) * ps[i]))
        adj[i] = run
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("res")
    ap.add_argument("ref")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--md")
    ap.add_argument("--title", default="")
    ap.add_argument("--holm", default="")
    ap.add_argument("--methods", default="")
    a = ap.parse_args()
    R, T, err = load(a.res, a.threads)
    methods = a.methods.split(",") if a.methods else sorted(R, key=lambda m: sum(R[m].values()) / len(R[m]))
    out = []
    P = out.append
    P(f"### {a.title}\n")
    P("**Mean FN rate per method** (n = datasets; sorted)\n")
    P("| method | n | all | FastMulRFS 25 bp | FastMulRFS 100 bp | DISCO | mean wall (s) | mean CPU (s) |")
    P("|---|---|---|---|---|---|---|---|")
    for m in methods:
        v = R[m]

        def mean(f):
            x = [v[k] for k in v if f(k)]
            return f"{sum(x) / len(x):.4f}" if x else "–"
        tw = [T[m][k][0] for k in v]
        tc = [T[m][k][1] for k in v]
        P(f"| {m} | {len(v)} | {mean(lambda k: True)} | {mean(lambda k: k[0] == 'fmrfs' and k[3] == 25)} | "
          f"{mean(lambda k: k[0] == 'fmrfs' and k[3] == 100)} | {mean(lambda k: k[0] == 'disco')} | "
          f"{sum(tw) / len(tw):.2f} | {sum(tc) / len(tc):.2f} |")
    P("")
    P(f"**Paired: {a.ref} minus method** (negative = {a.ref} better; W = {a.ref} better)\n")
    P("| vs | n | mean diff | 95% CI | W/T/L | p (Wilcoxon) | Holm p |")
    P("|---|---|---|---|---|---|---|")
    rows = [(m, compare(R, a.ref, m)) for m in methods if m != a.ref]
    rows = [(m, c) for m, c in rows if c]
    hm = [m for m in a.holm.split(",") if m]
    hp = dict(zip(hm, holm([dict(rows)[m]["p"] for m in hm]))) if hm and all(m in dict(rows) for m in hm) else {}
    for m, c in rows:
        h = f"{hp[m]:.2g}" if m in hp else ""
        P(f"| {m} | {c['n']} | {c['mean']:+.4f} | [{c['lo']:+.4f}, {c['hi']:+.4f}] | {c['W']}/{c['T']}/{c['L']} | {c['p']:.2g} | {h} |")
    P("")
    conds = sorted({(k[0], k[1], k[3]) for m in methods for k in R[m]})
    P("**Mean FN rate by condition**\n")
    P("| data | cond | sqln | " + " | ".join(methods) + " |")
    P("|---|---|---|" + "---|" * len(methods))
    for c in conds:
        cells = []
        vals = {}
        for m in methods:
            x = [R[m][k] for k in R[m] if (k[0], k[1], k[3]) == c]
            vals[m] = (sum(x) / len(x), len(x)) if x else None
        best = min(v[0] for v in vals.values() if v)
        for m in methods:
            v = vals[m]
            cells.append("–" if not v else (f"**{v[0]:.3f}**" if v[0] <= best + 1e-12 else f"{v[0]:.3f}") + (f" ({v[1]})" if v and v[1] != max(x[1] for x in vals.values() if x) else ""))
        P(f"| {c[0]} | {c[1]} | {c[2]} | " + " | ".join(cells) + " |")
    if err:
        P("\n**Failed runs**\n")
        for (m, d, c, e), n in sorted(err.items()):
            P(f"- {m} {d} {c}: {e} × {n}")
    s = "\n".join(out) + "\n"
    print(s)
    if a.md:
        with open(a.md, "a") as f:
            f.write(s + "\n")


if __name__ == "__main__":
    main()
