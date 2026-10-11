"""Convex 'node-measure' PDD (motivated by Theorem 1) vs DecoDiPhy on identical noisy inputs.

For each completed ksel run (same simulated d-hat, same pruned tree):
  stage 1:  RSS0 = min ||d - F m - ybar 1||^2  s.t. m >= 0, 1'm = 1, ybar >= 0   (m on ALL nodes)
  stage 2:  max ybar  s.t. same constraints and ||d - F m - ybar 1||^2 <= (1 + tau) RSS0
(two convex programs; no search over edge sets, no k). Accuracy = tree earth-mover distance (EMD)
between the true node measure and the estimate, EMD = sum_e l_e |M_true(below e) - M_hat(below e)|,
i.e. weighted UniFrac on placements with pendant lengths ignored (as in the paper's evaluation).
DecoDiPhy solutions (paper's k rule, and oracle k) are converted to node measures the same way
(a query on the edge above node a at relative position x from a puts p(1-x) on a and p x on parent(a)).
usage: python noisy_qp.py KSEL_OUT OUT.md
"""
import sys, os, glob, json, math, time
import numpy as np
import cvxpy as cp
from treeswift import read_tree_newick
from scipy.stats import wilcoxon

D, OUT = sys.argv[1], sys.argv[2]
MAXSEED = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
TRAIN = {"bees", "birds-jarvis", "1kp", "beetles", "hemipteroid"}
TAUS = [0.0, 0.05, 0.2, 0.5, 1.0]


def tree_mats(nwk, leaves_order=None):
    t = read_tree_newick(nwk)
    nodes = list(t.traverse_postorder())
    idx = {nd: i for i, nd in enumerate(nodes)}
    leaves = [nd for nd in nodes if nd.is_leaf()]
    V = len(nodes)
    # distances from every leaf to every node via BFS on the undirected tree
    adj = [[] for _ in range(V)]
    for nd in nodes:
        if nd.parent is not None:
            l = nd.edge_length or 0.0
            adj[idx[nd]].append((idx[nd.parent], l)); adj[idx[nd.parent]].append((idx[nd], l))
    F = np.zeros((len(leaves), V))
    for i, lf in enumerate(leaves):
        s = idx[lf]; dist = {s: 0.0}; st = [s]
        while st:
            u = st.pop()
            for v, l in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + l; st.append(v)
        F[i] = [dist[j] for j in range(V)]
    label = {nd.label: idx[nd] for nd in nodes}
    parent = [idx[nd.parent] if nd.parent is not None else -1 for nd in nodes]
    blen = [nd.edge_length or 0.0 for nd in nodes]
    return F, [lf.label for lf in leaves], label, parent, blen, V


def emd(m1, m2, parent, blen):
    diff = m1 - m2
    sub = diff.copy()
    tot = 0.0
    for i in range(len(sub)):  # postorder: children before parents
        if parent[i] >= 0:
            tot += blen[i] * abs(sub[i])
            sub[parent[i]] += sub[i]
    return tot


def to_measure(anchors, p, x, label, parent, V):
    m = np.zeros(V)
    for a, pq, xq in zip(anchors, p, x):
        i = label[a]
        m[i] += pq * (1 - xq)
        m[parent[i] if parent[i] >= 0 else i] += pq * xq
    return m


def solve(F, d, taus):
    n, V = F.shape
    m = cp.Variable(V, nonneg=True); y = cp.Variable(nonneg=True)
    r = d - F @ m - y
    cons = [cp.sum(m) == 1]
    sc = 1.0 / max(np.abs(d).max(), 1e-9)
    p1 = cp.Problem(cp.Minimize(cp.sum_squares(r * sc)), cons)
    p1.solve(solver="CLARABEL")
    rss0 = p1.value
    out = {}
    for tau in taus:
        p2 = cp.Problem(cp.Maximize(y), cons + [cp.sum_squares(r * sc) <= (1 + tau) * rss0 + 1e-12])
        try:
            p2.solve(solver="CLARABEL")
            out[tau] = (np.maximum(m.value, 0), float(y.value))
        except Exception:
            out[tau] = None
    return out


def rule_paper(r, thr=0.01):
    for rd in r["rounds"]:
        if min(rd["p"]) < thr:
            return max(rd["k"] - 1, 1)
        if rd["loss"] < 1e-10:
            return rd["k"]
    return r["rounds"][-1]["k"]


rows = []
cache = OUT + ".cache.jsonl"
done = {}
if os.path.exists(cache):
    for l in open(cache):
        z = json.loads(l); done[z["id"]] = z
for fn in sorted(glob.glob(os.path.join(D, "*.json"))):
    r = json.load(open(fn))
    if r["seed"] > MAXSEED:
        continue
    rid = os.path.basename(fn)
    if rid in done:
        rows.append(done[rid]); continue
    sim = fn[:-5] + "_sim"
    F, leaves, label, parent, blen, V = tree_mats(open(os.path.join(sim, "pruned_tree.trees")).read().strip())
    dist = {a: float(b) for a, b in (l.split() for l in open(os.path.join(sim, "distances.txt")))}
    d = np.array([dist[l] for l in leaves])
    true = [l.split() for l in open(os.path.join(sim, "true_queries.txt"))]
    mt = to_measure([t[0] for t in true], [float(t[1]) for t in true], [float(t[2]) for t in true], label, parent, V)
    t0 = time.time()
    sols = solve(F, d, TAUS)
    tq = time.time() - t0
    row = dict(id=rid, tree=r["tree"], k=r["k"], noise=r["noise"], seed=r["seed"], qp_time=tq,
               dd_time=sum(rd["time"] for rd in r["rounds"] if rd["k"] <= rule_paper(r) + 1))
    kp = rule_paper(r)
    for name, kk in [("decodiphy_paper", kp), ("decodiphy_oraclek", r["k"])]:
        rd = next(x for x in r["rounds"] if x["k"] == kk)
        row[name] = emd(mt, to_measure(rd["anchors"], rd["p"], rd["x"], label, parent, V), parent, blen)
    for tau, s in sols.items():
        row[f"qp_tau{tau}"] = emd(mt, s[0], parent, blen) if s is not None else float("nan")
        row[f"qp_tau{tau}_support"] = int((s[0] > 0.01).sum()) if s is not None else -1
    rows.append(row)
    with open(cache, "a") as f:
        f.write(json.dumps(row) + "\n")

train = [x for x in rows if x["tree"] in TRAIN]
best_tau = min(TAUS, key=lambda t: np.nanmean([x[f"qp_tau{t}"] for x in train]))
lines = [f"runs: {len(rows)}; tau tuned on train trees: {best_tau}",
         "", "| set | noise | n | DecoDiPhy (paper k) EMD | DecoDiPhy (oracle k) EMD | convex max-ybar EMD | dEMD (convex - paper) | Wilcoxon p | convex better / worse | median support (m>0.01) vs true 2k |",
         "|---|---|---|---|---|---|---|---|---|---|"]
for setname, S in [("test", [x for x in rows if x["tree"] not in TRAIN]), ("all", rows)]:
    for noise in ["noise0", "noise1", "noise2"]:
        sub = [x for x in S if x["noise"] == noise]
        if not sub:
            continue
        a = np.array([x["decodiphy_paper"] for x in sub]); o = np.array([x["decodiphy_oraclek"] for x in sub])
        q = np.array([x[f"qp_tau{best_tau}"] for x in sub])
        diff = q - a
        nz = diff[np.abs(diff) > 1e-9]
        p = wilcoxon(nz).pvalue if len(nz) >= 5 else float("nan")
        sup = np.median([x[f"qp_tau{best_tau}_support"] for x in sub]); tk = np.median([2 * x["k"] for x in sub])
        lines.append(f"| {setname} | {noise} | {len(sub)} | {a.mean():.4f} | {o.mean():.4f} | {q.mean():.4f} | {diff.mean():+.4f} | {p:.2g} | {(nz < 0).sum()}/{(nz > 0).sum()} | {sup:.0f} vs {tk:.0f} |")
lines.append("")
lines.append(f"mean runtime per instance: convex (2 + {len(TAUS)-1} solves, all taus) {np.mean([x['qp_time'] for x in rows]):.2f}s; "
             f"DecoDiPhy greedy up to chosen k+1 {np.mean([x['dd_time'] for x in rows]):.2f}s")
for t in TAUS:
    lines.append(f"tau={t}: mean EMD all={np.nanmean([x[f'qp_tau{t}'] for x in rows]):.4f}")
txt = "\n".join(lines)
print(txt)
open(OUT, "w").write(txt + "\n")
