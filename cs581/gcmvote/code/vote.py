# AI-assisted (Claude), exploration code for CS581 project
"""Vote model for MAGUS's GCM evidence: reference-free posterior that a cross-subset edge is true.

    python vote.py REPDIR VARIANT [--bb DIR] [--out OUT.fasta] [--tag NAME]   -> merged alignment

REPDIR is a replicate in the bbevidence/gcmgen layout: REPDIR/inputs/subalignments (MAGUS's 25 subset
alignments) and REPDIR/inputs/backbones (MAGUS's 10 L-INS-i backbones, backbone_N_mafft.txt). --bb uses
another backbone directory (e.g. 5 or 20 backbones). Output: REPDIR/vote/<tag>/out.fasta (or --out),
plus edges.npz (every cross-subset edge: nodes, k, n, weight, posterior) and model.json (fit, cutoff).

For every cross-subset GCM edge (column a of subset A, column b of subset B):
  k = support  = number of backbones that align at least one residue pair between a and b;
  n = exposure = number of backbones that COULD vote: they hold a residue of A in column a and a
                 residue of B in column b (the backbone covers both nodes).
Only edges with k >= 1 are observed, so the model is a zero-truncated two-component mixture,
fitted by EM per run, without the reference:
  true edges  k | n ~ Bin(n, p1),   false edges  k | n ~ Bin(n, p0),   P(true) = pi    (binomial)
  or the same with BetaBin(n, a, b) components                                       (beta-binomial)
Posterior w = P(true | k, n).

VARIANT:
  magus            MAGUS's own graph (control; must reproduce MAGUS's merge)
  esK              delete edges with support < K (gcmgen's hand-picked filter is es4)
  fracF            delete edges with support < ceil(F * B), B = number of backbones (e.g. frac0.4)
  hard[-bb]        keep edges with w > 0.5             (binomial / beta-binomial mixture)
  soft[-bb]        edge weight * w
  softG[-bb]       edge weight * w^G  (e.g. soft2, soft4-bb)
  hardT[-bb]       keep edges with w > T (e.g. hard0.9)
  ...+mask         tree variant: additionally write out.masked.fasta with the final columns whose
                   cross-subset evidence is mostly low-posterior removed (see mask_columns)
Weights are kept integer (x1000, rounded) because MAGUS's graph file is integer; MCL is invariant to a
global scaling of the weights (checked by gcmgen with file duplication).
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
from scipy.optimize import minimize, minimize_scalar
from scipy.special import betaln, gammaln

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, "..", "..", "code"))
MERGE_FLAGS = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
               "-f", "4"]


# ---------------------------------------------------------------- mixture models on (k, n) tables

def _lchoose(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def _binom_logpmf_trunc(k, n, p):
    p = min(max(p, 1e-9), 1 - 1e-9)
    return _lchoose(n, k) + k * np.log(p) + (n - k) * np.log1p(-p) - np.log1p(-(1 - p) ** n)


def _bb_logpmf_trunc(k, n, a, b):
    lp = _lchoose(n, k) + betaln(k + a, n - k + b) - betaln(a, b)
    p0 = np.exp(betaln(a, n + b) - betaln(a, b))  # P(k = 0)
    return lp - np.log1p(-np.minimum(p0, 1 - 1e-12))


def fit_mixture(K, N, C, kind="binom", iters=500, tol=1e-9):
    """EM on the aggregated table: cells (K[i], N[i]) observed C[i] times (K >= 1).
    Returns dict(pi, params of both components, loglik, bic) and the per-cell posterior P(true)."""
    K, N, C = (np.asarray(x, dtype=float) for x in (K, N, C))
    tot = C.sum()
    if kind == "binom":
        th1, th0 = [0.8], [0.15]
        logpmf = lambda th: _binom_logpmf_trunc(K, N, th[0])
    else:  # beta-binomial: (mean, log concentration)
        th1, th0 = [0.8, math.log(5.0)], [0.15, math.log(5.0)]

        def logpmf(th):
            m, s = th[0], math.exp(th[1])
            return _bb_logpmf_trunc(K, N, m * s, (1 - m) * s)
    pi, ll_old = 0.5, -np.inf
    for it in range(iters):
        l1 = np.log(pi) + logpmf(th1)
        l0 = np.log(1 - pi) + logpmf(th0)
        mx = np.maximum(l1, l0)
        lse = mx + np.log(np.exp(l1 - mx) + np.exp(l0 - mx))
        r = np.exp(l1 - lse)
        ll = float((C * lse).sum())
        pi = float(np.clip((C * r).sum() / tot, 1e-6, 1 - 1e-6))
        th1 = _mstep(K, N, C * r, th1, kind)
        th0 = _mstep(K, N, C * (1 - r), th0, kind)
        if abs(ll - ll_old) < tol * abs(ll):
            break
        ll_old = ll
    if th1[0] < th0[0]:  # label switching: component 1 = the high-vote (true) component
        th1, th0, pi, r = th0, th1, 1 - pi, 1 - r
    npar = 3 if kind == "binom" else 5
    res = {"kind": kind, "pi": pi, "p1": th1[0], "p0": th0[0], "loglik": ll, "iters": it + 1,
           "bic": -2 * ll + npar * math.log(tot), "n_edges": int(tot)}
    if kind == "bb":
        res.update(conc1=math.exp(th1[1]), conc0=math.exp(th0[1]))
    return res, r


def _mstep(K, N, Wt, th, kind):
    if Wt.sum() <= 0:
        return th
    if kind == "binom":
        f = lambda p: -(Wt * _binom_logpmf_trunc(K, N, p)).sum()
        return [minimize_scalar(f, bounds=(1e-6, 1 - 1e-6), method="bounded").x]

    def f(x):
        m, s = 1 / (1 + math.exp(-x[0])), math.exp(x[1])
        return -(Wt * _bb_logpmf_trunc(K, N, m * s, (1 - m) * s)).sum()
    x0 = [math.log(th[0] / (1 - th[0])), th[1]]
    x = minimize(f, x0, method="L-BFGS-B", bounds=[(-12, 12), (-5, 8)]).x
    return [1 / (1 + math.exp(-x[0])), x[1]]


def table(k, n):
    """Aggregate edges into (k, n) cells. Returns K, N, C and the cell index of every edge."""
    key = n.astype(np.int64) * 1000 + k
    uk, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    return uk % 1000, uk // 1000, cnt, inv


def fit_gof(K, N, C, fit):
    """Expected counts of every (k, n) cell under the fitted mixture (given the observed n margin),
    and the G statistic / df."""
    exp = np.zeros(len(C))
    for nn in np.unique(N):
        sel = N == nn
        ks = np.arange(1, nn + 1)
        if fit["kind"] == "binom":
            p = fit["pi"] * np.exp(_binom_logpmf_trunc(ks, nn, fit["p1"])) + \
                (1 - fit["pi"]) * np.exp(_binom_logpmf_trunc(ks, nn, fit["p0"]))
        else:
            s1, s0 = fit["conc1"], fit["conc0"]
            p = fit["pi"] * np.exp(_bb_logpmf_trunc(ks, nn, fit["p1"] * s1, (1 - fit["p1"]) * s1)) + \
                (1 - fit["pi"]) * np.exp(_bb_logpmf_trunc(ks, nn, fit["p0"] * s0, (1 - fit["p0"]) * s0))
        # (the mixture's n-specific mixing weight is not modelled; renormalise within n)
        p = p / p.sum()
        exp[sel] = C[sel].sum() * p[(K[sel] - 1).astype(int)]
    ok = C > 0
    G = 2 * float((C[ok] * np.log(C[ok] / np.maximum(exp[ok], 1e-12))).sum())
    return exp, G


# ---------------------------------------------------------------- graph construction (MAGUS-exact)

def build(context):
    """MAGUS's graph (sum over backbones of A_b^T A_b, as gcmx.fastgraph) plus support and exposure."""
    from magus.configuration import Configs
    from magus.tasks import task
    from gcmx.fastgraph import _column_matrix
    from gcmx.weighting import _backbone_alignmap
    graph = context.graph
    n = graph.matrixSize
    files = [t.outputFile for t in task.asCompleted(context.backboneTasks)]
    for path in context.backbonePaths:
        if path not in files:
            files.append(path)
    files = sorted(files)
    total = sp.csr_matrix((n, n), dtype=np.int64)
    support = sp.csr_matrix((n, n), dtype=np.int64)
    cover = np.zeros(n, dtype=np.int64)  # bit b set: backbone b holds a residue in node
    for bi, f in enumerate(files):
        A = _column_matrix(_backbone_alignmap(context, f), n)
        prod = (A.T @ A).tocsr()
        total = total + prod
        support = support + (prod > 0).astype(np.int64)
        cover |= (np.asarray(A.sum(axis=0)).ravel() > 0).astype(np.int64) << bi
    Configs.log("[vote] {} backbones".format(len(files)))
    return total.tocsr(), support.tocsr(), cover, len(files)


def popcount(x):
    x = x.astype(np.uint64)
    c = np.zeros(len(x), dtype=np.int64)
    while x.any():
        c += (x & np.uint64(1)).astype(np.int64)
        x >>= np.uint64(1)
    return c


def edge_arrays(context, total, support, cover):
    """Cross-subset edges with a < b: rows, cols, weight, k, n."""
    sub = np.array([s for s, _ in context.graph.matSubPosMap])
    coo = sp.triu(total, k=1).tocoo()
    keep = sub[coo.row] != sub[coo.col]
    r, c, w = coo.row[keep], coo.col[keep], coo.data[keep]
    k = np.asarray(support[r, c]).ravel()
    nexp = popcount(cover[r] & cover[c])
    return r, c, w, k, nexp, sub


def transform(context, total, support, cover, nbb, variant, outdir):
    from magus.configuration import Configs
    r, c, w, k, nexp, sub = edge_arrays(context, total, support, cover)
    assert (k >= 1).all() and (nexp >= k).all(), "exposure must bound support"
    base = variant.replace("+mask", "")
    model = {"variant": variant, "B": nbb, "edges": int(len(k))}
    post = np.ones(len(k))
    kind = "bb" if base.endswith("-bb") else "binom"
    K, N, C, inv = table(k, nexp)
    fits = {}
    for kd in ("binom", "bb"):
        fit, rc = fit_mixture(K, N, C, kd)
        exp, G = fit_gof(K, N, C, fit)
        fit["G"], fit["cells"] = G, int(len(C))
        fits[kd] = (fit, rc)
    fit, rcell = fits[kind]
    post = rcell[inv]
    model["fit"] = {kd: f for kd, (f, _) in fits.items()}
    # effective cutoff: smallest k with posterior > 0.5, per exposure n
    model["cutoff_by_n"] = {}
    for nn in range(1, nbb + 1):
        sel = N == nn
        if sel.any():
            ks = K[sel][rcell[sel] > 0.5]
            model["cutoff_by_n"][nn] = int(ks.min()) if len(ks) else None
    model["n_hist"] = {int(nn): int(C[N == nn].sum()) for nn in np.unique(N)}
    np.savez_compressed(os.path.join(outdir, "edges.npz"), a=r, b=c, sub_a=sub[r], sub_b=sub[c], w=w, k=k, n=nexp,
                        post=post, post_binom=fits["binom"][1][inv], post_bb=fits["bb"][1][inv])
    # new weights
    name = base.replace("-bb", "")
    if name == "magus":
        neww = w.astype(float)
    elif name.startswith("es"):
        neww = np.where(k >= int(name[2:]), w, 0)
    elif name.startswith("frac"):
        neww = np.where(k >= math.ceil(float(name[4:]) * nbb - 1e-9), w, 0)
    elif name.startswith("hard"):
        t = float(name[4:]) if len(name) > 4 else 0.5
        neww = np.where(post > t, w, 0)
    elif name.startswith("soft"):
        g = float(name[4:]) if len(name) > 4 else 1.0
        neww = w * post ** g
    else:
        raise SystemExit("unknown variant " + variant)
    model["kept_edges"] = int((neww > 0).sum())
    model["kept_weight_frac"] = float(neww.sum() / w.sum())
    json.dump(model, open(os.path.join(outdir, "model.json"), "w"), indent=1, default=float)
    if name == "magus":
        return total
    scale = 1 if name.startswith(("es", "frac", "hard")) else 1000
    neww = np.rint(neww * scale).astype(np.int64)
    # rebuild: non-cross entries (self loops, same-subset) keep weight * scale; cross edges get neww (symmetric)
    coo = total.tocoo()
    cross = sub[coo.row] != sub[coo.col]
    keep_r, keep_c, keep_v = coo.row[~cross], coo.col[~cross], coo.data[~cross] * scale
    nz = neww > 0
    rows = np.concatenate([keep_r, r[nz], c[nz]])
    cols = np.concatenate([keep_c, c[nz], r[nz]])
    vals = np.concatenate([keep_v, neww[nz], neww[nz]])
    Configs.log("[vote] {}: kept {} of {} cross edges".format(variant, int(nz.sum()), len(nz)))
    return sp.csr_matrix((vals, (rows, cols)), shape=total.shape)


def install(variant, outdir):
    import numpy as _np
    from magus.align.merge.graph_build import graph_builder as gb
    from magus.configuration import Configs
    from gcmx import fastgraph

    def buildMatrix(context):
        total, support, cover, nbb = build(context)
        total = transform(context, total, support, cover, nbb, variant, outdir).tocsr()
        if Configs.graphBuildRestrict:
            raise SystemExit("graphbuildrestrict not supported")
        total.indices = total.indices.astype(_np.int32)
        context.graph.matrix = fastgraph.CSRGraph(total)
    fastgraph.install()
    gb.buildMatrix = buildMatrix


# ---------------------------------------------------------------- column mask for trees

def mask_columns(outdir, aln_path, sub_dir, thresh=0.5):
    """Remove final columns whose cross-subset evidence is mostly low-posterior.
    For a final column, the evidence is every cross-subset edge (from the full, unfiltered graph) whose two
    nodes both landed in that column; score = weight-weighted mean posterior of those edges. Columns
    with score < thresh are removed; columns with no cross-subset edge (single-subset columns) are kept."""
    sys.path.insert(0, CODE)
    from gcmx import fasta
    E = np.load(os.path.join(outdir, "edges.npz"))
    files = sorted(os.listdir(sub_dir), key=lambda f: int(f.split("_")[-1].split(".")[0]) if f.split("_")[-1].split(".")[0].isdigit() else f)
    aln = fasta.read(aln_path)
    # node id = subsetMatrixIdx[sub] + subset column. MAGUS orders subsets as context.subalignmentPaths
    order = json.load(open(os.path.join(outdir, "subsets.json")))
    off, nodecol = {}, {}
    o = 0
    node_final = []
    for si, path in enumerate(order):
        sa = fasta.read(path)
        L = len(next(iter(sa.values())))
        col = np.full(L, -1, dtype=np.int64)
        for t, s in sa.items():  # residues of t: subset column -> final column
            sc = [i for i, ch in enumerate(s) if ch not in "-."]
            fc = [i for i, ch in enumerate(aln[t]) if ch not in "-."]
            col[np.array(sc, dtype=np.int64)] = np.array(fc, dtype=np.int64)
        node_final.append(col)
    node_final = np.concatenate(node_final)
    fa, fb = node_final[E["a"]], node_final[E["b"]]
    same = (fa == fb) & (fa >= 0)
    L = len(next(iter(aln.values())))
    wsum = np.bincount(fa[same], weights=E["w"][same], minlength=L)
    psum = np.bincount(fa[same], weights=(E["w"] * E["post"])[same], minlength=L)
    score = np.where(wsum > 0, psum / np.maximum(wsum, 1e-12), 1.0)
    keep = score >= thresh
    out = {t: "".join(ch for ch, kk in zip(s, keep) if kk) for t, s in aln.items()}
    dst = aln_path.replace(".fasta", ".masked.fasta")
    fasta.write(out, dst)
    return dst, int((~keep).sum()), L


# ---------------------------------------------------------------- driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rep")
    ap.add_argument("variant")
    ap.add_argument("--bb", help="backbone directory (default REP/inputs/backbones)")
    ap.add_argument("--out")
    ap.add_argument("--tag")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--keep-work", action="store_true")
    a = ap.parse_args()
    rep = os.path.abspath(a.rep)
    tag = a.tag or a.variant
    vd = os.path.join(rep, "vote", tag)
    shutil.rmtree(vd, ignore_errors=True)
    os.makedirs(vd)
    bb = os.path.abspath(a.bb or os.path.join(rep, "inputs", "backbones"))
    out = os.path.abspath(a.out) if a.out else os.path.join(vd, "out.fasta")
    subdir = os.path.join(rep, "inputs", "subalignments")
    start = time.time()
    with open(os.path.join(vd, "magus.log"), "w") as log:
        subprocess.run([sys.executable, os.path.abspath(__file__), "--_inner", a.variant, vd,
                        "-np", str(a.threads), "-d", os.path.join(vd, "work"), "-s", subdir, "-b", bb, "-o", out]
                       + MERGE_FLAGS, cwd=CODE, stdout=log, stderr=subprocess.STDOUT, check=True)
    wall = round(time.time() - start, 1)
    if not a.keep_work:
        shutil.rmtree(os.path.join(vd, "work"), ignore_errors=True)
    info = {"merge_wall": wall, "out": out}
    if a.variant.endswith("+mask"):
        dst, nmask, L = mask_columns(vd, out, subdir)
        info.update(masked=dst, masked_cols=nmask, cols=L)
    json.dump(info, open(os.path.join(vd, "run.json"), "w"))
    print(json.dumps(info))
    return out


def inner():
    variant, vd = sys.argv[2], sys.argv[3]
    sys.argv = [sys.argv[0]] + sys.argv[4:]
    sys.path.insert(0, CODE)
    install(variant, vd)
    # record MAGUS's subset order (node ids are offsets in this order)
    from magus.align.merge.graph_build import graph_builder as gb
    orig = gb.buildMatrix

    def wrapped(context):
        json.dump([os.path.abspath(p) for p in context.subalignmentPaths], open(os.path.join(vd, "subsets.json"), "w"))
        return orig(context)
    gb.buildMatrix = wrapped
    from magus.main import main as magus_main
    magus_main()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--_inner":
        inner()
    else:
        main()
