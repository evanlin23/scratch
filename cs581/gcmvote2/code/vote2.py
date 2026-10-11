# AI-assisted (Claude), exploration code for CS581 project
"""Vote models v2 for MAGUS's GCM evidence: graded, reliability-weighted and overlap-aware votes.

    python vote2.py REPDIR VARIANT [--bb DIR] [--tag NAME]   -> REPDIR/vote2/<tag>/out.fasta

Same harness as gcmvote/code/vote.py (vote_v1.py here): MAGUS's own merge (-f 4, MCL, minclusters) on its
own subsets and backbones; only the cross-subset edge weights of the alignment graph change.
The first run on a replicate caches the per-backbone evidence in REPDIR/vote2/features[_<bbtag>].npz;
later variants reuse it (the graph is rebuilt bit-for-bit from the cache).

Per cross-subset edge e = (column a of subset A, column b of subset B) and backbone j:
  x_j  = residue pairs backbone j aligns between a and b (MAGUS's weight is w = sum_j x_j)
  r_j(a), r_j(b) = residues backbone j holds in node a / b
  exposed_j = r_j(a) > 0 and r_j(b) > 0;   n = #exposed;   k = #{j: x_j > 0}
  s_j  = x_j / (r_j(a) r_j(b))   graded vote: the share of the residue pairs backbone j could align
         between a and b that it does align (in [0, 1]; 0 for an exposed non-vote).

VARIANT (hard = keep edges with posterior > 0.5; "-soft" = weight x posterior):
  magus, esK, hard-bb     baselines (hard-bb = gcmvote v1: zero-truncated beta-binomial mixture on (k, n))
  m1-gbb                  M1: beta-binomial mixture on the graded support k_s = round(sum_j s_j) in [1, n]
  m1-dm[-soft]            M1: Dirichlet-multinomial mixture on (#none, #weak (s <= .5), #strong (s > .5))
  m2-ds[-soft]            M2: binary Dawid-Skene (per-backbone sensitivity / false-vote rate), zero-truncated
  m2-dsg                  M2: graded Dawid-Skene (per backbone, categorical none / weak / strong)
  m3-ovbb                 M3: hard-bb fit, per-edge log-likelihood ratio x (n_eff/n) / median(n_eff/n) (n_eff from
                          the overlap of the backbones' sequence sets for subsets A, B);
                          guard: if the edge-median n_eff < 2, the votes are not independent -> magus
  m23-ovds                M2+M3: Dawid-Skene log-likelihood ratio tempered by n_eff / n; same guard
  m4-gmm[-soft]           M4: 2-component full-covariance Gaussian mixture on
                          (logit (k+.5)/(n+1), logit mean s over voting backbones, log w)
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
from scipy.special import gammaln, expit, logit

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, "..", "..", "code"))
sys.path.insert(0, HERE)
import vote_v1 as v1  # noqa: E402

MERGE_FLAGS = v1.MERGE_FLAGS
EPS = 1e-6


# ---------------------------------------------------------------- feature extraction (inside MAGUS)

def extract(context, path):
    """Union graph (MAGUS-exact) + per-backbone pairs / residues on every cross-subset edge."""
    from magus.tasks import task
    from gcmx.fastgraph import _column_matrix
    from gcmx.weighting import _backbone_alignmap
    from magus.helpers import sequenceutils
    graph = context.graph
    n = graph.matrixSize
    files = [t.outputFile for t in task.asCompleted(context.backboneTasks)]
    for p in context.backbonePaths:
        if p not in files:
            files.append(p)
    files = sorted(files)
    sub = np.array([s for s, _ in graph.matSubPosMap])
    mats, res, taxa = [], [], []
    total = sp.csr_matrix((n, n), dtype=np.int64)
    for f in files:
        A = _column_matrix(_backbone_alignmap(context, f), n)
        prod = (A.T @ A).tocsr()
        total = total + prod
        mats.append(prod)
        res.append(np.asarray(A.sum(axis=0)).ravel())
        taxa.append(list(sequenceutils.readFromFasta(f).keys()))
    total = total.tocsr()
    coo = sp.triu(total, k=1).tocoo()
    keep = sub[coo.row] != sub[coo.col]
    r, c, w = coo.row[keep], coo.col[keep], coo.data[keep]
    key = r.astype(np.int64) * n + c
    order = np.argsort(key)
    r, c, w, key = r[order], c[order], w[order], key[order]
    B = len(files)
    X = np.zeros((len(r), B), dtype=np.int32)
    for j, prod in enumerate(mats):
        pc = sp.triu(prod, k=1).tocoo()
        kk = pc.row.astype(np.int64) * n + pc.col
        pos = np.searchsorted(key, kk)
        ok = (pos < len(key)) & (key[np.minimum(pos, len(key) - 1)] == kk)
        X[pos[ok], j] = pc.data[ok]
    RA = np.stack([rr[r] for rr in res], axis=1).astype(np.int32)
    RB = np.stack([rr[c] for rr in res], axis=1).astype(np.int32)
    # backbone sequence sets per subset (for overlap / n_eff)
    tsub = {t: s for t, s in context.taxonSubalignmentMap.items()}
    nsub = int(sub.max()) + 1
    members = [[set() for _ in range(nsub)] for _ in range(B)]
    for j, tl in enumerate(taxa):
        for t in tl:
            if t in tsub:
                members[j][tsub[t]].add(t)
    nsets = np.array([[len(members[j][s]) for s in range(nsub)] for j in range(B)])
    inter = np.zeros((B, B, nsub), dtype=np.int64)
    for i in range(B):
        for j in range(B):
            for s in range(nsub):
                inter[i, j, s] = len(members[i][s] & members[j][s])
    np.savez_compressed(path, r=r, c=c, w=w, sub_a=sub[r], sub_b=sub[c], X=X, RA=RA, RB=RB, nsets=nsets,
                        inter=inter, nodesub=sub, t_data=total.data, t_indices=total.indices, t_indptr=total.indptr,
                        n=n, files=np.array(files))


def load(path):
    F = dict(np.load(path, allow_pickle=False))
    F["total"] = sp.csr_matrix((F["t_data"], F["t_indices"], F["t_indptr"]), shape=(int(F["n"]), int(F["n"])))
    E = F["RA"] > 0
    E &= F["RB"] > 0
    F["E"] = E
    F["V"] = F["X"] > 0
    F["k"] = F["V"].sum(1)
    F["nexp"] = E.sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        S = F["X"] / np.maximum(F["RA"].astype(float) * F["RB"], 1)
    F["S"] = np.where(E, np.clip(S, 0, 1), 0.0)
    return F


# ---------------------------------------------------------------- models

def m_hardbb(F):
    return _bb_fit_int(F["k"], F["nexp"])


def _bb_fit_int(k, n):
    """vote_v1's zero-truncated beta-binomial mixture on integer (k, n) with 1 <= k <= n."""
    K, N, C, inv = v1.table(k.astype(np.int64), n.astype(np.int64))
    fit, rc = v1.fit_mixture(K, N, C, "bb")
    return rc[inv], fit


def m1_gbb(F):
    """Graded support k_s = sum_j s_j, rounded to an integer in [1, n] (a continuous-k beta-binomial is not a
    proper likelihood: it degenerated on the first fit check)."""
    n = F["nexp"]
    ks = np.clip(np.rint(F["S"].sum(1)), 1, n)
    return _bb_fit_int(ks, n)


def _dm_logpmf(M, a):
    """log DM pmf of count rows M (cells x 3) with concentration a (3,), and log P(all in category 0)."""
    n = M.sum(1)
    A = a.sum()
    lp = gammaln(n + 1) - gammaln(M + 1).sum(1) + gammaln(A) - gammaln(n + A) + (gammaln(M + a) - gammaln(a)).sum(1)
    lp0 = gammaln(A) - gammaln(n + A) + gammaln(n + a[0]) - gammaln(a[0])
    return lp, lp0


def m1_dm(F, iters=300, tol=1e-9):
    E, S = F["E"], F["S"]
    m1 = ((S > 0) & (S <= 0.5)).sum(1)
    m2 = (S > 0.5).sum(1)
    m0 = F["nexp"] - m1 - m2
    key = (m0 * 100 + m1) * 100 + m2
    uk, inv, C = np.unique(key, return_inverse=True, return_counts=True)
    M = np.stack([uk // 10000, (uk // 100) % 100, uk % 100], axis=1).astype(float)
    C = C.astype(float)

    def ll_trunc(a):
        lp, lp0 = _dm_logpmf(M, a)
        return lp - np.log1p(-np.minimum(np.exp(lp0), 1 - 1e-12))
    th = [np.log([1.0, 1.0, 4.0]), np.log([4.0, 1.0, 0.5])]
    pi, old = 0.5, -np.inf
    for it in range(iters):
        l1 = math.log(pi) + ll_trunc(np.exp(th[0]))
        l0 = math.log(1 - pi) + ll_trunc(np.exp(th[1]))
        mx = np.maximum(l1, l0)
        lse = mx + np.log(np.exp(l1 - mx) + np.exp(l0 - mx))
        r = np.exp(l1 - lse)
        ll = float((C * lse).sum())
        pi = float(np.clip((C * r).sum() / C.sum(), 1e-6, 1 - 1e-6))
        for ci, wt in ((0, C * r), (1, C * (1 - r))):
            f = lambda x, wt=wt: -(wt * ll_trunc(np.exp(x))).sum()
            th[ci] = minimize(f, th[ci], method="L-BFGS-B", bounds=[(-6, 8)] * 3).x
        if abs(ll - old) < tol * abs(ll):
            break
        old = ll
    a1, a0 = np.exp(th[0]), np.exp(th[1])
    if a1[2] / a1.sum() < a0[2] / a0.sum():
        a1, a0, pi, r = a0, a1, 1 - pi, 1 - r
    return r[inv], {"pi": pi, "alpha_true": a1.tolist(), "alpha_false": a0.tolist(), "loglik": ll, "iters": it + 1}


def _patterns(F, graded):
    """Unique per-backbone observation patterns: state 0 unexposed, 1 exposed no vote, 2 vote (binary) or
    2 weak / 3 strong (graded). Returns pattern matrix (P x B), counts, inverse."""
    E, S = F["E"], F["S"]
    st = np.where(E, 1, 0)
    if graded:
        st = np.where(E & (S > 0) & (S <= 0.5), 2, st)
        st = np.where(E & (S > 0.5), 3, st)
    else:
        st = np.where(E & F["V"], 2, st)
    base = 4 if graded else 3
    code = (st.astype(np.int64) * base ** np.arange(st.shape[1])).sum(1)
    uc, inv, cnt = np.unique(code, return_inverse=True, return_counts=True)
    P = np.stack([(uc // base ** j) % base for j in range(st.shape[1])], axis=1)
    return P, cnt.astype(float), inv


def ds_fit(F, graded=False, iters=500, tol=1e-10):
    """Dawid-Skene latent class model with per-backbone class-conditional vote distributions.
    Edges with no vote are unobserved (zero truncation): EM with the unobserved all-'no vote' edges as missing
    data (expected count per class and exposure mask), so pi is the class share among all exposed node pairs."""
    P, C, inv = _patterns(F, graded)
    B = P.shape[1]
    ncat = 3 if graded else 2  # categories among exposed: 0 = no vote, 1.. = vote levels
    Ex = P > 0
    Y = np.zeros(P.shape + (ncat,))
    for cat in range(ncat):
        Y[..., cat] = (P == cat + 1)
    # init: true votes often, false rarely
    th = np.zeros((2, B, ncat))
    if graded:
        th[0] = [0.15, 0.25, 0.6]
        th[1] = [0.8, 0.15, 0.05]
    else:
        th[0] = [0.2, 0.8]
        th[1] = [0.8, 0.2]
    pi = np.array([0.1, 0.9])
    old = -np.inf
    for it in range(iters):
        lt = np.log(np.clip(th, 1e-6, 1))
        L = np.einsum("pbk,cbk->pc", Y, lt)          # log P(observed pattern | class)
        lp0 = np.einsum("pb,cb->pc", Ex.astype(float), lt[:, :, 0])  # log P(no vote at all | class, mask)
        lj = np.log(pi)[None, :] + L
        mx = lj.max(1, keepdims=True)
        lse = mx[:, 0] + np.log(np.exp(lj - mx).sum(1))
        R = np.exp(lj - lse[:, None])
        # observed-data loglik: sum C log( sum_c pi_c P(y|c) / sum_c pi_c (1 - P0_c) )
        pobs = (pi[None, :] * (1 - np.exp(lp0))).sum(1)
        ll = float((C * (lse - np.log(np.maximum(pobs, 1e-300)))).sum())
        Rc = C[:, None] * R
        p0 = np.exp(lp0)
        U = Rc * p0 / np.maximum(1 - p0, 1e-12)      # expected unobserved edges with this mask, per class
        for ci in range(2):
            num = np.einsum("p,pbk->bk", Rc[:, ci], Y)
            num[:, 0] += np.einsum("p,pb->b", U[:, ci], Ex.astype(float))
            den = num.sum(1, keepdims=True)
            th[ci] = np.clip(num / np.maximum(den, 1e-12), 1e-4, 1 - 1e-4)
            th[ci] /= th[ci].sum(1, keepdims=True)
        tot = Rc.sum(0) + U.sum(0)
        pi = np.clip(tot / tot.sum(), 1e-9, 1 - 1e-9)
        if abs(ll - old) < tol * abs(ll):
            break
        old = ll
    # true class = the one with the higher mean vote probability
    vote1 = (1 - th[:, :, 0]).mean(1)
    tcls = int(np.argmax(vote1))
    lt = np.log(np.clip(th, 1e-6, 1))
    L = np.einsum("pbk,cbk->pc", Y, lt)
    llr = (L[:, tcls] - L[:, 1 - tcls])
    prior = math.log(pi[tcls] / pi[1 - tcls])
    post = expit(prior + llr)
    fit = {"pi_true_all_exposed": float(pi[tcls]), "iters": it + 1, "loglik": ll, "graded": graded,
           "true_cat_probs": th[tcls].tolist(), "false_cat_probs": th[1 - tcls].tolist(),
           "post_obs_share": float((C * post).sum() / C.sum())}
    return post[inv], fit, (prior, llr[inv])


def m2_ds(F):
    p, fit, _ = ds_fit(F, graded=False)
    return p, fit


def m2_dsg(F):
    p, fit, _ = ds_fit(F, graded=True)
    return p, fit


def neff(F):
    """Effective number of independent votes per edge: |S|^2 / sum_{i,j in S} J_ij, where J_ij is the Jaccard
    overlap of the sequences of subsets A and B held by backbones i and j (S = exposed backbones)."""
    inter, nsets = F["inter"], F["nsets"]
    sa, sb, E = F["sub_a"], F["sub_b"], F["E"]
    B = E.shape[1]
    mask = (E.astype(np.int64) << np.arange(B)).sum(1)
    pairkey = (sa.astype(np.int64) * 10000 + sb) * (1 << B) + mask
    uk, inv = np.unique(pairkey, return_inverse=True)
    out = np.zeros(len(uk))
    for u, kk in enumerate(uk):
        m = int(kk % (1 << B))
        ab = int(kk // (1 << B))
        a, b = ab // 10000, ab % 10000
        idx = [j for j in range(B) if m >> j & 1]
        I = inter[np.ix_(idx, idx)][:, :, a] + inter[np.ix_(idx, idx)][:, :, b]
        sz = nsets[idx, a] + nsets[idx, b]
        U = sz[:, None] + sz[None, :] - I
        J = np.where(U > 0, I / np.maximum(U, 1), 0.0)
        np.fill_diagonal(J, 1.0)
        out[u] = len(idx) ** 2 / J.sum()
    return out[inv]


def m3_ovbb(F):
    """hard-bb fit on (k, n); per-edge log-likelihood ratio tempered by t = (n_eff/n) / median(n_eff/n), clipped
    to [0, 2]: the beta-binomial already absorbs the average vote correlation, so only edges whose backbones
    overlap more (less) than typical are discounted (boosted)."""
    ne = neff(F)
    med = float(np.median(ne))
    if med < 2:
        return None, {"guard": True, "median_neff": med}
    p, fit = _bb_fit_int(F["k"], F["nexp"])
    pi = fit["pi"]
    lr = logit(np.clip(p, 1e-12, 1 - 1e-12)) - math.log(pi / (1 - pi))
    rho = ne / F["nexp"]
    t = np.clip(rho / np.median(rho), 0, 2)
    fit.update(guard=False, median_neff=med, t_quantiles=np.quantile(t, [0.05, 0.25, 0.5, 0.75, 0.95]).tolist())
    return expit(math.log(pi / (1 - pi)) + t * lr), fit


def m23_ovds(F):
    ne = neff(F)
    med = float(np.median(ne))
    if med < 2:
        return None, {"guard": True, "median_neff": med}
    _, fit, (prior, llr) = ds_fit(F, graded=False)
    t = ne / F["nexp"]
    fit.update(guard=False, median_neff=med)
    return expit(prior + t * llr), fit


def m4_gmm(F, nfit=300000, seed=0):
    from sklearn.mixture import GaussianMixture
    k, n, S, V = F["k"], F["nexp"], F["S"], F["V"]
    ms = (S * V).sum(1) / np.maximum(k, 1)
    Z = np.stack([logit((k + 0.5) / (n + 1.0)), logit(np.clip(ms, 0.01, 0.99)), np.log(F["w"].astype(float))], 1)
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(Z), min(nfit, len(Z)), replace=False)
    g = GaussianMixture(2, covariance_type="full", n_init=3, random_state=seed).fit(Z[idx])
    tc = int(np.argmax(g.means_[:, 0]))
    p = g.predict_proba(Z)[:, tc]
    return p, {"means": g.means_.tolist(), "weights": g.weights_.tolist(), "true_comp": tc}


MODELS = {"hard-bb": m_hardbb, "m1-gbb": m1_gbb, "m1-dm": m1_dm, "m2-ds": m2_ds, "m2-dsg": m2_dsg,
          "m3-ovbb": m3_ovbb, "m23-ovds": m23_ovds, "m4-gmm": m4_gmm}


def weights(F, variant):
    """New cross-edge weights (float) and model info."""
    w, k = F["w"].astype(float), F["k"]
    if variant == "magus":
        return w, {}
    if variant.startswith("es"):
        return np.where(k >= int(variant[2:]), w, 0.0), {}
    soft = variant.endswith("-soft")
    base = variant[:-5] if soft else variant
    post, info = MODELS[base](F)
    if post is None:  # guard: fall back to MAGUS
        return w, info
    info["kept_edges"] = int((post > 0.5).sum())
    info["post_hist"] = np.histogram(post, bins=10, range=(0, 1))[0].tolist()
    return (w * post if soft else np.where(post > 0.5, w, 0.0)), info


def rebuild(F, neww, soft):
    total = F["total"]
    sub = F["nodesub"]
    scale = 1000 if soft else 1
    nw = np.rint(neww * scale).astype(np.int64)
    coo = total.tocoo()
    cross = sub[coo.row] != sub[coo.col]
    kr, kc, kv = coo.row[~cross], coo.col[~cross], coo.data[~cross] * scale
    nz = nw > 0
    rows = np.concatenate([kr, F["r"][nz], F["c"][nz]])
    cols = np.concatenate([kc, F["c"][nz], F["r"][nz]])
    vals = np.concatenate([kv, nw[nz], nw[nz]])
    return sp.csr_matrix((vals, (rows, cols)), shape=total.shape)


def install(variant, vd, cache):
    import numpy as _np
    from magus.align.merge.graph_build import graph_builder as gb
    from magus.configuration import Configs
    from gcmx import fastgraph

    def buildMatrix(context):
        json.dump([os.path.abspath(p) for p in context.subalignmentPaths], open(os.path.join(vd, "subsets.json"), "w"))
        if not os.path.exists(cache):
            extract(context, cache + ".tmp.npz")
            os.rename(cache + ".tmp.npz", cache)
        F = load(cache)
        assert F["total"].shape[0] == context.graph.matrixSize
        t0 = time.time()
        neww, info = weights(F, variant)
        info.update(variant=variant, B=int(F["X"].shape[1]), edges=int(len(neww)),
                    kept_edges_final=int((neww > 0).sum()), kept_weight_frac=float(neww.sum() / F["w"].sum()),
                    fit_seconds=round(time.time() - t0, 1))
        json.dump(info, open(os.path.join(vd, "model.json"), "w"), indent=1, default=float)
        np.save(os.path.join(vd, "neww.npy"), neww.astype(np.float32))
        total = F["total"] if variant == "magus" else rebuild(F, neww, variant.endswith("-soft"))
        if Configs.graphBuildRestrict:
            raise SystemExit("graphbuildrestrict not supported")
        total = total.tocsr()
        total.indices = total.indices.astype(_np.int32)
        context.graph.matrix = fastgraph.CSRGraph(total)
        Configs.log("[vote2] {}: kept {} of {} cross edges".format(variant, info["kept_edges_final"], len(neww)))
    fastgraph.install()
    gb.buildMatrix = buildMatrix


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rep")
    ap.add_argument("variant")
    ap.add_argument("--bb")
    ap.add_argument("--tag")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    rep = os.path.abspath(a.rep)
    tag = a.tag or a.variant
    vd = os.path.join(rep, "vote2", tag)
    shutil.rmtree(vd, ignore_errors=True)
    os.makedirs(vd)
    bb = os.path.abspath(a.bb or os.path.join(rep, "inputs", "backbones"))
    cache = os.path.join(rep, "vote2", "features_{}.npz".format(os.path.basename(bb)))
    out = os.path.join(vd, "out.fasta")
    subdir = os.path.join(rep, "inputs", "subalignments")
    start = time.time()
    with open(os.path.join(vd, "magus.log"), "w") as log:
        subprocess.run([sys.executable, os.path.abspath(__file__), "--_inner", a.variant, vd, cache,
                        "-np", str(a.threads), "-d", os.path.join(vd, "work"), "-s", subdir, "-b", bb, "-o", out]
                       + MERGE_FLAGS, cwd=CODE, stdout=log, stderr=subprocess.STDOUT, check=True)
    shutil.rmtree(os.path.join(vd, "work"), ignore_errors=True)
    info = {"merge_wall": round(time.time() - start, 1), "out": out}
    json.dump(info, open(os.path.join(vd, "run.json"), "w"))
    print(json.dumps(info))


def inner():
    variant, vd, cache = sys.argv[2], sys.argv[3], sys.argv[4]
    sys.argv = [sys.argv[0]] + sys.argv[5:]
    sys.path.insert(0, CODE)
    install(variant, vd, cache)
    from magus.main import main as magus_main
    magus_main()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--_inner":
        inner()
    else:
        main()
