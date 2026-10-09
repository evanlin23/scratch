"""Estimate the (haploid) effective population size Ne of the Nute et al. simulations from
rooted-triplet concordance of true gene trees: P(gene triplet = species triplet) = 1 - 2/3 exp(-t/Ne)."""
import sys, os, numpy as np, treeswift as ts
from scipy.optimize import minimize_scalar
ROOT = "/opt/data/nute"
obs = []  # (t_generations, n_match, n_total)
for h in ["10M", "2M", "500K"]:
    for r in ["1E-7", "1E-6"]:
        for rep in range(1, 6):
            d = f"{ROOT}/25tax-1000gen-0bps-{h}-{r}-full/{rep:02d}"
            sp = ts.read_tree_newick(open(f"{d}/true-species.tre").read().strip())
            trip = []
            for u in sp.traverse_internal():
                if u.parent is None or len(u.children) != 2: continue
                sib = [c for c in u.parent.children if c is not u]
                a = next(u.children[0].traverse_leaves()).label
                b = next(u.children[1].traverse_leaves()).label
                c = next(sib[0].traverse_leaves()).label
                trip.append((a, b, c, u.edge_length))
            genes = [ts.read_tree_newick(l) for l in open(f"{d}/true-genes.tre").read().split("\n")[:300] if l.strip()]
            for a, b, c, t in trip:
                m = 0
                for g in genes:
                    lab = {l.label: l for l in g.traverse_leaves()}
                    anc = set()
                    x = lab[a]
                    while x is not None:
                        anc.add(id(x)); x = x.parent
                    # LCA(a,b) depth vs LCA(a,c)
                    def lca(y):
                        y = lab[y]
                        while id(y) not in anc: y = y.parent
                        return y
                    lab_ab, lab_ac = lca(b), lca(c)
                    # ab|c iff LCA(a,b) is strictly below LCA(a,c)
                    y = lab_ab; below = False
                    while y is not None:
                        if y is lab_ac and lab_ab is not lab_ac: below = True; break
                        y = y.parent
                    m += below
                obs.append((h, t, m, len(genes)))
obs_t = np.array([o[1] for o in obs]); M = np.array([o[2] for o in obs]); N = np.array([o[3] for o in obs])
def nll(logNe):
    P = np.clip(1 - 2/3 * np.exp(-obs_t / np.exp(logNe)), 1e-9, 1 - 1e-9)
    return -(M * np.log(P) + (N - M) * np.log(1 - P)).sum()
res = minimize_scalar(nll, bounds=(np.log(1e3), np.log(1e8)), method="bounded")
Ne = np.exp(res.x)
print("Ne_hat =", round(Ne), " (#edges", len(obs), ")")
for h in ["10M", "2M", "500K"]:
    sel = np.array([o[0] == h for o in obs])
    r2 = minimize_scalar(lambda l: -((M[sel]*np.log(np.clip(1-2/3*np.exp(-obs_t[sel]/np.exp(l)),1e-9,1)))+(N[sel]-M[sel])*np.log(np.clip(2/3*np.exp(-obs_t[sel]/np.exp(l)),1e-9,1))).sum(), bounds=(np.log(1e3), np.log(1e8)), method="bounded")
    print(h, "Ne_hat", round(np.exp(r2.x)))
