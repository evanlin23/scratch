"""Exact interior edge lengths beta(c) of the limiting ASTRID-Pro metric (theory.md, Thm 3-4).

(1) The SimPhy GDL settings of the DISCO and FastMulRFS benchmarks (constant per-generation rates on
    the true species trees): is the theorem's condition (all interior beta > 0) met?
(2) Random per-branch rates on random Yule trees with 4-8 taxa: how often is some interior beta < 0,
    and are all such cases supercritical interior branches above lossy subtrees?
Usage: python beta_scan.py OUT.json"""
import glob
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gdl", "code"))
sys.path.insert(0, HERE)
import gdlsim  # noqa: E402
from phylo import parse_newick  # noqa: E402
from verify_theory import betas, survival, mean_copies  # noqa: E402


def interior_betas(st, s):
    b = betas(st, s)
    return {c: v for c, v in b.items() if st.parent[st.parent[c]] >= 0}


def bench_settings():
    out = {}
    disco = {"default": (5e-10, 1), "gdl_1e-10_1": (1e-10, 1), "gdl_1e-10_05": (1e-10, .5),
             "gdl_1e-10_0": (1e-10, 0), "gdl_5e-10_05": (5e-10, .5), "gdl_5e-10_0": (5e-10, 0),
             "gdl_1e-9_1": (1e-9, 1), "gdl_1e-9_05": (1e-9, .5), "gdl_1e-9_0": (1e-9, 0)}
    for c, (lam, ratio) in disco.items():
        mins = []
        for f in sorted(glob.glob(f"/opt/data/disco/trees/{c}/*/s_tree.trees")):
            st = parse_newick(open(f).read().strip())
            n = len(st.parent)
            s = survival(st, [lam] * n, [lam * ratio] * n)
            mins.append(round(min(interior_betas(st, s).values()), 4))
        out["disco/" + c] = {"lambda": lam, "mu": lam * ratio, "min_interior_beta_per_rep": mins}
    for dl in ("0.0000000001", "0.0000000002", "0.0000000005"):
        mins = []
        for f in sorted(glob.glob(f"/opt/data/fmrfs/ntaxa-100.dlrate-{dl}.psize-10000000/*/s_tree.trees")):
            st = parse_newick(open(f).read().strip())
            n = len(st.parent)
            s = survival(st, [float(dl)] * n, [float(dl)] * n)
            mins.append(round(min(interior_betas(st, s).values()), 4))
        out["fmrfs/dl" + dl] = {"lambda": float(dl), "mu": float(dl), "min_interior_beta_per_rep": mins}
    return out


def random_scan(ntrials=20000, seed=1):
    rng = random.Random(seed)
    neg = 0
    neg_supercrit = 0
    neg_examples = []
    by_n = {}
    for t in range(ntrials):
        k = rng.randint(4, 8)
        st = parse_newick(gdlsim.yule_tree(k, rng.randrange(10**9), height=rng.choice([1.0, 2.0, 3.0])))
        n = len(st.parent)
        lam = [rng.choice([0, 0.5, 1, 2, 3]) for _ in range(n)]
        mu = [rng.choice([0, 0.5, 1, 2, 3]) for _ in range(n)]
        s = survival(st, lam, mu)
        ib = interior_betas(st, s)
        by_n.setdefault(k, [0, 0])
        by_n[k][1] += 1
        if ib and min(ib.values()) < 0:
            neg += 1
            by_n[k][0] += 1
            c = min(ib, key=ib.get)
            sc = mean_copies(lam[c], mu[c], st.length[c]) > 1
            neg_supercrit += sc
            if len(neg_examples) < 5:
                neg_examples.append({"tree": st.newick(lengths=True), "c": c, "beta": round(ib[c], 4),
                                     "lam_c": lam[c], "mu_c": mu[c]})
    return {"trials": ntrials, "neg": neg, "neg_with_supercritical_branch": neg_supercrit,
            "by_ntaxa": {k: {"neg": v[0], "trials": v[1]} for k, v in sorted(by_n.items())},
            "examples": neg_examples}


if __name__ == "__main__":
    res = {"benchmarks": bench_settings(), "random": random_scan()}
    json.dump(res, open(sys.argv[1], "w"), indent=1)
    for k, v in res["benchmarks"].items():
        print(k, "min beta over reps:", min(v["min_interior_beta_per_rep"]))
    r = res["random"]
    print("random: %d/%d negative, %d of them with supercritical branch" % (r["neg"], r["trials"],
                                                                          r["neg_with_supercritical_branch"]))
    print(r["by_ntaxa"])
