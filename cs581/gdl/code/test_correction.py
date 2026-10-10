"""Test ASTRID-Pro-S (survival-reweighted) on configurations where ideal ASTRID-Pro fails.
First pass: (a) oracle = true rooted species tree; (b) the plain ASTRID-Pro tree (true tags,
root counted) rooted at the given outgroup, then iterate (re-estimate s_hat on the new tree).
Usage: python test_correction.py TREE RATES_JSON NFAM SEED OUTGROUP [min_species]"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import methods as M  # noqa: E402
from correction import SpTree, corrected_distances, estimate_survival, reroot_species_tree  # noqa: E402
from gdlsim import simulate  # noqa: E402
from phylo import parse_newick, rf_error  # noqa: E402

tree, rates_s, nfam, seed, outg = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
minsp = int(sys.argv[6]) if len(sys.argv) > 6 else 2
st = parse_newick(tree)
rates = json.loads(rates_s)
lam = [rates.get(st.label[v], [0, 0])[0] for v in range(len(st.parent))]
mu = [rates.get(st.label[v], [0, 0])[1] for v in range(len(st.parent))]
species = sorted(st.label[v] for v in st.leaves())
idx = {s: i for i, s in enumerate(species)}
fams, _, _ = simulate(st, lam, mu, nfam, seed, min_species=minsp)
G = [parse_newick(f) for f in fams]
tagged = [M.root_and_tag(g, M.simphy_species, idx, truetags=True) for g in G]
res = {"tree": tree, "nfam": nfam, "seed": seed}
per = [M.gene_distances(rt, tg, ls, len(species), mode="pro", count_root=True) for rt, tg, ls in tagged]
D, _ = M.average_matrix(per)
T0 = M.fastme_tree(D, species)
res["pro_true_root_FN"] = rf_error(T0, st)[0]
per = [M.gene_distances(rt, tg, ls, len(species), mode="multi") for rt, tg, ls in tagged]
Tm = M.fastme_tree(M.average_matrix(per)[0], species)
res["multi_FN"] = rf_error(Tm, st)[0]


def corrected_tree(S):
    sh = estimate_survival(tagged, S)
    per = [corrected_distances(rt, tg, ls, len(species), S, sh) for rt, tg, ls in tagged]
    D, _ = M.average_matrix(per)
    return M.fastme_tree(D, species), sh


S_or = SpTree(st, idx)
T, sh = corrected_tree(S_or)
res["proS_oracle_FN"] = rf_error(T, st)[0]
res["shat_oracle"] = {(st.label[c] or "".join(sorted(st.label[x] for x in st.leaves() if (S_or.mask[c] >> idx[st.label[x]]) & 1))): round(v, 4)
                      for c, v in sh.items()}
cur = T0
hist = []
for it in range(4):
    S = SpTree(reroot_species_tree(cur, outg), idx)
    cur, _ = corrected_tree(S)
    hist.append(rf_error(cur, st)[0])
res["proS_iter_FN"] = hist
S = SpTree(reroot_species_tree(Tm, outg), idx)
T, _ = corrected_tree(S)
res["proS_from_multi_FN"] = rf_error(T, st)[0]
print(json.dumps(res), flush=True)
