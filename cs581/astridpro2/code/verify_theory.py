"""Check the additivity theorem (theory.md) numerically.

For a species tree with per-branch linear birth-death rates, compute the exact survival
probabilities s(c) (Kendall's generating function), the predicted limiting ASTRID-Pro matrix
delta(A,B) = 1 + sum_{v strictly inside path, v != LCA} s(off(v)), and the predicted edge lengths
beta(c) = (s(c1) + s(c2) + s(sib c) - s(c)) / 2. Then simulate families with true root and tags,
run apro -M pro -T (root counted) and compare the empirical matrix with delta.
Usage: python verify_theory.py OUT.jsonl"""
import json
import math
import os
import random
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gdl", "code"))
import gdlsim  # noqa: E402
from phylo import parse_newick  # noqa: E402

APRO = os.path.join(HERE, "apro")


def pgf(lam, mu, t, z):
    """E[z^N_t] for linear birth-death started from one copy."""
    if t <= 0 or (lam == 0 and mu == 0):
        return z
    if abs(lam - mu) < 1e-12:
        a = lam * t / (1 + lam * t)
        b = a
    else:
        e = math.exp((lam - mu) * t)
        a = mu * (e - 1) / (lam * e - mu)
        b = lam * (e - 1) / (lam * e - mu)
    return a + (1 - a) * (1 - b) * z / (1 - b * z)


def mean_copies(lam, mu, t):
    return math.exp((lam - mu) * t)


def survival(st, lam, mu):
    """s[v]: a copy entering the branch above v leaves >= 1 descendant at the leaves."""
    s = [0.0] * len(st.parent)
    for v in st.postorder():
        if st.children[v]:
            q = 1.0
            for c in st.children[v]:
                q *= 1 - s[c]
            q = 1 - q
        else:
            q = 1.0
        t = st.length[v] or 0.0
        s[v] = 1 - pgf(lam[v], mu[v], t, 1 - q)
    return s


def predicted(st, s, count_root=True):
    leaves = st.leaves()
    names = [st.label[v] for v in leaves]
    anc = {}
    for v in leaves:
        p, path = v, []
        while p >= 0:
            path.append(p)
            p = st.parent[p]
        anc[v] = path
    D = np.zeros((len(leaves), len(leaves)))
    for i, a in enumerate(leaves):
        for j, b in enumerate(leaves):
            if i >= j:
                continue
            sa = set(anc[b])
            L = next(x for x in anc[a] if x in sa)
            tot = 1.0
            for leaf in (a, b):
                prev = leaf
                for x in anc[leaf][1:]:
                    if x == L:
                        break
                    off = [c for c in st.children[x] if c != prev][0]
                    tot += s[off]
                    prev = x
            D[i, j] = D[j, i] = tot
    return names, D


def betas(st, s):
    out = {}
    for c in range(len(st.parent)):
        p = st.parent[c]
        if p < 0 or not st.children[c]:
            continue
        sib = [x for x in st.children[p] if x != c][0]
        c1, c2 = st.children[c]
        out[c] = 0.5 * (s[c1] + s[c2] + s[sib] - s[c])
    return out


FASTME = "/opt/mm/root/envs/gdl/bin/fastme"


def fn_of(D, names, st):
    """FastME (BalME + NNI + SPR) on matrix D; FN against the species tree st."""
    from phylo import rf_error
    with tempfile.TemporaryDirectory() as td:
        i, o = os.path.join(td, "d.phy"), os.path.join(td, "t.nwk")
        with open(i, "w") as f:
            f.write("%d\n" % len(names))
            for a, nm in enumerate(names):
                f.write(nm + " " + " ".join("%.6f" % x for x in D[a]) + "\n")
        subprocess.run([FASTME, "-i", i, "-o", o, "-m", "B", "-n", "B", "-s", "-T", "1"], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return rf_error(parse_newick(open(o).read()), st)[0]


def empirical(fams, names, args=("-M", "pro", "-T", "-R", "1"), first=None):
    with tempfile.TemporaryDirectory() as td:
        g, o = os.path.join(td, "g.nwk"), os.path.join(td, "o.phy")
        open(g, "w").write("\n".join(fams) + "\n")
        extra = []
        if first:
            open(os.path.join(td, "s.nwk"), "w").write(first + "\n")
            extra = ["-s", os.path.join(td, "s.nwk")]
        subprocess.run([APRO, "-i", g, "-o", o, "-u"] + list(args) + extra, check=True,
                       stderr=subprocess.DEVNULL)
        L = [l.split() for l in open(o).read().split("\n")[1:] if l.strip()]
    idx = {r[0]: i for i, r in enumerate(L)}
    D = np.full((len(names), len(names)), np.nan)
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if a in idx and b in idx:
                D[i, j] = float(L[idx[a]][1 + idx[b]])
    return D


def run(cfg, nfam, seed, out):
    st = parse_newick(cfg["tree"])
    n = len(st.parent)
    lam, mu = [0.0] * n, [0.0] * n
    for v in range(n):
        key = st.label[v] if not st.children[v] else None
        rates = cfg["rates"].get(str(v)) or (cfg["rates"].get(key) if key else None) or cfg["default"]
        lam[v], mu[v] = rates
    s = survival(st, lam, mu)
    names, P = predicted(st, s)
    fams, tries, over = gdlsim.simulate(st, lam, mu, nfam, seed, min_species=2, cap=4000)
    E = empirical(fams, names)
    iu = np.triu_indices(len(names), 1)
    err = np.nanmax(np.abs(E[iu] - P[iu]))
    b = betas(st, s)
    fns = {"pred": fn_of(P, names, st), "pro_truetags": fn_of(E, names, st)}
    for nm, args in (("pro_inferred", ("-M", "pro")), ("multi", ("-M", "multi"))):
        fns[nm] = fn_of(empirical(fams, names, args), names, st)
    # Pro-S with true tags; first pass = the true (rooted) species tree, and = the unrooted multi tree
    strip = lambda t: t  # noqa: E731
    fns["proS_truetags_oracle"] = fn_of(empirical(fams, names, ("-M", "pros", "-T"), first=cfg["tree"]), names, st)
    fns["proS_inferred_oracle"] = fn_of(empirical(fams, names, ("-M", "pros"), first=cfg["tree"]), names, st)
    rec = {"name": cfg["name"], "nfam": nfam, "seed": seed, "max_abs_err": round(float(err), 4),
           "min_beta": round(min(b.values()), 4), "FN": fns, "betas": {str(k): round(v, 4) for k, v in b.items()},
           "s": [round(x, 4) for x in s], "overflow": over,
           "pred": {f"{names[i]}-{names[j]}": round(float(P[i, j]), 4) for i, j in zip(*iu)},
           "emp": {f"{names[i]}-{names[j]}": round(float(E[i, j]), 4) for i, j in zip(*iu)}}
    print(json.dumps({k: rec[k] for k in ("name", "nfam", "seed", "max_abs_err", "min_beta", "FN")}), flush=True)
    with open(out, "a") as f:
        f.write(json.dumps(rec) + "\n")


def label_tree(nwk):
    return nwk


if __name__ == "__main__":
    out = sys.argv[1]
    # node ids follow phylo.parse_newick preorder; rates keyed by leaf label or node id
    cfgs = [
        {"name": "counterexample", "tree": "(((A:1,B:1):1,C:1):1,D:3);",
         "rates": {"A": [0, 3], "B": [0, 3], "C": [0, 3], "2": [3, 0]}, "default": [0, 0]},
        {"name": "critical_6taxon", "tree": "(((A:0.5,B:0.5):0.4,C:0.9):0.3,((D:0.6,E:0.6):0.3,F:0.9):0.3);",
         "rates": {}, "default": [1.0, 1.0]},
        {"name": "subcritical_lossy", "tree": "(((A:0.5,B:0.5):0.4,C:0.9):0.3,((D:0.6,E:0.6):0.3,F:0.9):0.3);",
         "rates": {}, "default": [0.5, 1.5]},
    ]
    # random supercritical configurations on a 6-taxon caterpillar
    rng = random.Random(1)
    for r in range(4):
        rates = {}
        for v in range(11):
            rates[str(v)] = [rng.choice([0, 0.5, 1, 2, 3]), rng.choice([0, 0.5, 1, 2, 3])]
        cfgs.append({"name": f"random{r}", "tree": "(((((A:0.7,B:0.7):0.4,C:1.1):0.3,D:1.4):0.5,E:1.9):0.2,F:2.1);",
                     "rates": rates, "default": [0, 0]})
    nfam = int(os.environ.get("NFAM", "20000"))
    for c in cfgs:
        run(c, nfam, 1, out)
