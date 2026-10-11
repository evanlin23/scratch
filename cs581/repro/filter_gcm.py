# AI-assisted (Claude), code for CS581 project
"""Filter low-support evidence in MAGUS's merge step (GCM), then let MAGUS finish the merge.

Merge only (the subsets and backbones already exist, e.g. a bank replicate):

    python3 filter_gcm.py WORKDIR VARIANT --out merged.fasta [--ref true.fasta] [--backbones DIR] [-np 4]

    WORKDIR/inputs/subalignments/subalignment_subset_*.txt   MAGUS's subset alignments
    WORKDIR/inputs/backbones/backbone_*_mafft.txt            MAGUS's backbone alignments
    WORKDIR/subsets.json  (optional)                         subset order of the logged run (see subset_order)

End to end (MAGUS decomposes, aligns subsets and backbones itself, then merges with the filter switched on):

    python3 filter_gcm.py --from-scratch UNALIGNED.fasta VARIANT --out merged.fasta [--subsets 25] [--ref ...]

VARIANT (all of them keep MAGUS's MCL clustering and min-clusters trace unchanged; only the graph changes):
    magus           MAGUS's own graph (control)
    esK             drop every cross-subset edge supported by fewer than K backbones (es4 = K 4)
    fracF           drop edges supported by fewer than ceil(F * B) backbones, B = number of backbones
    vote-hard       keep edges whose posterior P(true | k, n) > 0.5, binomial mixture        (= vote.py "hard")
    vote-hard-bb    the same with a beta-binomial mixture (the pre-registered selected variant, = "hard-bb")
    vote-hardT[-bb] keep edges with posterior > T (e.g. vote-hard0.9)

With --ref the merged alignment is scored with FastSP (as gcmx.bbtool_bench.acc_ref does) and one JSON row is
printed: SPFN, SPFP, avgErr = (SPFN + SPFP) / 2, TC, LenEst, LenRef, kept_edges, edges, merge_wall.

Requirements: MAGUS (commit 39041fc, `pip install -e MAGUS`), numpy, scipy; for --ref: java and FastSP.jar
(path in $FASTSP_JAR, default /opt/tools/FastSP/FastSP.jar).

This is a cleaned-up copy of cs581/gcmvote/code/vote.py (branch claude/cs581-gcmvote) restricted to the variants
above; the graph and the edge weights it hands to MCL are byte-identical to vote.py's.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize, minimize_scalar
from scipy.special import betaln, gammaln

# MAGUS options for the merge, exactly as the logged runs used them (gcmx.bbtool_bench.MERGE_FLAGS)
MERGE_FLAGS = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
               "-f", "4"]


def scratch_flags(num_subsets):
    """MAGUS end-to-end options of the MAGUS paper (gcmx.e2e_bench.magus_flags): 25 subsets, 10 L-INS-i backbones
    of 200 sequences, MCL with inflation 4, min-clusters trace without the optimizer."""
    return ["--maxsubsetsize", "0", "--maxnumsubsets", str(num_subsets), "--decompstrategy", "pastastyle",
            "--decompskeletonsize", "300", "--graphbuildmethod", "mafft", "--graphclustermethod", "mcl",
            "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false", "-r", "10", "-m", "200", "-f", "4"]


# =====================================================================================================
# 1. The evidence: MAGUS's graph, plus how many backbones support / could support every edge
# =====================================================================================================

def backbone_files(context):
    """All backbone alignments MAGUS uses (the ones it computed itself plus the ones given with -b), sorted."""
    from magus.tasks import task
    files = [t.outputFile for t in task.asCompleted(context.backboneTasks)]
    files += [p for p in context.backbonePaths if p not in files]
    return sorted(files)


def column_matrix(context, backbone_path):
    """Sparse matrix A (backbone columns x graph nodes): A[l, v] = how many residues of node v sit in backbone
    column l. A graph node is one column of one subset alignment; MAGUS's own backboneToAlignMap finds it."""
    from magus.align.merge.graph_build import graph_builder as gb
    from magus.helpers import sequenceutils
    if backbone_path in context.backboneExtend:
        raise SystemExit("HMM-extended backbones are not supported (run MAGUS with --graphbuildhmmextend false)")
    aln = sequenceutils.readFromFasta(backbone_path)
    length = len(next(iter(aln.values())).seq)
    alignmap = gb.backboneToAlignMap(context, aln, length)  # list over columns of {node: count}
    rows, cols, vals = [], [], []
    for col, nodes in enumerate(alignmap):
        for node, count in nodes.items():
            rows.append(col)
            cols.append(node)
            vals.append(count)
    n = context.graph.matrixSize
    return sp.csr_matrix((np.array(vals, dtype=np.int64), (rows, cols)), shape=(len(alignmap), n))


def build_evidence(context):
    """Return (weight, support, cover, B).
    weight[a, b]  = MAGUS's edge weight: residue pairs between nodes a and b summed over backbones (= sum_b A^T A)
    support[a, b] = k, the number of backbones with at least one residue pair between a and b
    cover[v]      = bit mask of the backbones that hold at least one residue in node v
    B             = number of backbones"""
    n = context.graph.matrixSize
    weight = sp.csr_matrix((n, n), dtype=np.int64)
    support = sp.csr_matrix((n, n), dtype=np.int64)
    cover = np.zeros(n, dtype=np.int64)
    files = backbone_files(context)
    for i, path in enumerate(files):
        A = column_matrix(context, path)
        pairs = (A.T @ A).tocsr()                 # this backbone's residue pairs between every two nodes
        weight = weight + pairs
        support = support + (pairs > 0).astype(np.int64)
        holds = np.asarray(A.sum(axis=0)).ravel() > 0
        cover |= holds.astype(np.int64) << i
    return weight.tocsr(), support.tocsr(), cover, len(files)


def count_bits(x):
    """Number of set bits of every entry of an integer array."""
    x = x.astype(np.uint64)
    c = np.zeros(len(x), dtype=np.int64)
    while x.any():
        c += (x & np.uint64(1)).astype(np.int64)
        x >>= np.uint64(1)
    return c


def cross_edges(context, weight, support, cover):
    """Every edge between two different subsets, listed once (a < b), with
    w = weight, k = support, n = exposure = number of backbones holding a residue in BOTH nodes (the
    backbones that could have voted for the edge). Only edges with k >= 1 exist in the graph."""
    subset_of = np.array([s for s, _ in context.graph.matSubPosMap])
    upper = sp.triu(weight, k=1).tocoo()
    cross = subset_of[upper.row] != subset_of[upper.col]
    a, b, w = upper.row[cross], upper.col[cross], upper.data[cross]
    k = np.asarray(support[a, b]).ravel()
    n = count_bits(cover[a] & cover[b])
    assert (k >= 1).all() and (n >= k).all()
    return a, b, w, k, n, subset_of


# =====================================================================================================
# 2. The vote model: is an edge true, given that k of its n possible backbones voted for it?
# =====================================================================================================
# Two kinds of edges are mixed in the graph: true homologies (most backbones that could vote do vote) and false
# ones (few do). A two-component mixture is fitted to the (k, n) counts by EM, without any reference alignment:
#     true:  k | n ~ Bin(n, p1)      false:  k | n ~ Bin(n, p0)      P(true) = pi          ("binom")
# or the same with beta-binomial components BetaBin(n, mean, concentration)                ("bb").
# Edges with k = 0 never enter the graph, so both components are truncated at k >= 1.

def log_choose(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def binom_logpmf(k, n, p):
    """log P(k | n, p) for a binomial truncated to k >= 1."""
    p = min(max(p, 1e-9), 1 - 1e-9)
    return log_choose(n, k) + k * np.log(p) + (n - k) * np.log1p(-p) - np.log1p(-(1 - p) ** n)


def betabinom_logpmf(k, n, a, b):
    """log P(k | n, a, b) for a beta-binomial truncated to k >= 1."""
    lp = log_choose(n, k) + betaln(k + a, n - k + b) - betaln(a, b)
    p_zero = np.exp(betaln(a, n + b) - betaln(a, b))
    return lp - np.log1p(-np.minimum(p_zero, 1 - 1e-12))


def component_logpmf(K, N, theta, kind):
    """theta = [p] for binom, [mean, log concentration] for bb."""
    if kind == "binom":
        return binom_logpmf(K, N, theta[0])
    m, s = theta[0], math.exp(theta[1])
    return betabinom_logpmf(K, N, m * s, (1 - m) * s)


def refit_component(K, N, weights, theta, kind):
    """M step: maximise the weighted log-likelihood of one component."""
    if weights.sum() <= 0:
        return theta
    if kind == "binom":
        f = lambda p: -(weights * binom_logpmf(K, N, p)).sum()
        return [minimize_scalar(f, bounds=(1e-6, 1 - 1e-6), method="bounded").x]

    def f(x):  # x = (logit mean, log concentration)
        m, s = 1 / (1 + math.exp(-x[0])), math.exp(x[1])
        return -(weights * betabinom_logpmf(K, N, m * s, (1 - m) * s)).sum()
    x0 = [math.log(theta[0] / (1 - theta[0])), theta[1]]
    x = minimize(f, x0, method="L-BFGS-B", bounds=[(-12, 12), (-5, 8)]).x
    return [1 / (1 + math.exp(-x[0])), x[1]]


def fit_mixture(K, N, C, kind, iters=500, tol=1e-9):
    """EM on the table of distinct (k, n) cells; cell i was seen C[i] times.
    Returns the fitted parameters and, per cell, the posterior P(true | k, n)."""
    K, N, C = (np.asarray(x, dtype=float) for x in (K, N, C))
    theta1, theta0 = ([0.8], [0.15]) if kind == "binom" else ([0.8, math.log(5.0)], [0.15, math.log(5.0)])
    pi, ll_old = 0.5, -np.inf
    for it in range(iters):
        # E step: posterior responsibility of the "true" component for every cell
        l1 = np.log(pi) + component_logpmf(K, N, theta1, kind)
        l0 = np.log(1 - pi) + component_logpmf(K, N, theta0, kind)
        mx = np.maximum(l1, l0)
        lse = mx + np.log(np.exp(l1 - mx) + np.exp(l0 - mx))
        post = np.exp(l1 - lse)
        ll = float((C * lse).sum())
        # M step
        pi = float(np.clip((C * post).sum() / C.sum(), 1e-6, 1 - 1e-6))
        theta1 = refit_component(K, N, C * post, theta1, kind)
        theta0 = refit_component(K, N, C * (1 - post), theta0, kind)
        if abs(ll - ll_old) < tol * abs(ll):
            break
        ll_old = ll
    if theta1[0] < theta0[0]:  # name the components so that "true" is the one with more votes
        theta1, theta0, pi, post = theta0, theta1, 1 - pi, 1 - post
    fit = {"kind": kind, "pi": pi, "p1": theta1[0], "p0": theta0[0], "loglik": ll, "iters": it + 1}
    if kind == "bb":
        fit.update(conc1=math.exp(theta1[1]), conc0=math.exp(theta0[1]))
    return fit, post


def edge_posteriors(k, n, kind):
    """Fit the mixture on the (k, n) table of all cross-subset edges; return each edge's posterior and the fit."""
    key = n.astype(np.int64) * 1000 + k
    cells, cell_of_edge, count = np.unique(key, return_inverse=True, return_counts=True)
    K, N = cells % 1000, cells // 1000
    fit, post_cell = fit_mixture(K, N, count, kind)
    # effective cutoff: the smallest k the model keeps, for each exposure n
    fit["cutoff_by_n"] = {int(nn): (int(K[(N == nn) & (post_cell > 0.5)].min())
                                    if ((N == nn) & (post_cell > 0.5)).any() else None) for nn in np.unique(N)}
    return post_cell[cell_of_edge], fit


# =====================================================================================================
# 3. The variants: which cross-subset edges survive
# =====================================================================================================

def parse_variant(variant):
    """'vote-hard-bb' -> ('hard', 0.5, 'bb'); 'es4' -> ('es', 4, None); 'frac0.4' -> ('frac', 0.4, None)."""
    v = variant[5:] if variant.startswith("vote-") else variant
    kind = "bb" if v.endswith("-bb") else "binom"
    v = v[:-3] if v.endswith("-bb") else v
    if v == "magus":
        return "magus", None, None
    if v.startswith("es"):
        return "es", int(v[2:]), None
    if v.startswith("frac"):
        return "frac", float(v[4:]), None
    if v.startswith("hard"):
        return "hard", float(v[4:]) if len(v) > 4 else 0.5, kind
    raise SystemExit("unknown variant " + variant)


def keep_mask(variant, k, n, B, info):
    """Boolean array: which cross-subset edges keep their weight (the others are deleted)."""
    rule, x, kind = parse_variant(variant)
    if rule == "magus":
        return np.ones(len(k), dtype=bool)
    if rule == "es":
        return k >= x
    if rule == "frac":
        return k >= math.ceil(x * B - 1e-9)
    post, fit = edge_posteriors(k, n, kind)
    info["fit"] = fit
    return post > x


def filtered_graph(context, weight, support, cover, B, variant, info):
    """MAGUS's graph with the dropped cross-subset edges removed. Self loops and edges inside one subset are
    untouched; kept edges keep MAGUS's integer weight (so MCL sees exactly MAGUS's numbers on them)."""
    a, b, w, k, n, subset_of = cross_edges(context, weight, support, cover)
    keep = keep_mask(variant, k, n, B, info)
    info.update(B=B, edges=int(len(k)), kept_edges=int(keep.sum()),
                kept_weight_frac=float(w[keep].sum() / w.sum()))
    if parse_variant(variant)[0] == "magus":
        return weight
    coo = weight.tocoo()
    inside = subset_of[coo.row] == subset_of[coo.col]
    rows = np.concatenate([coo.row[inside], a[keep], b[keep]])
    cols = np.concatenate([coo.col[inside], b[keep], a[keep]])
    vals = np.concatenate([coo.data[inside], w[keep], w[keep]]).astype(np.int64)
    return sp.csr_matrix((vals, (rows, cols)), shape=weight.shape)


# =====================================================================================================
# 4. Plugging the filtered graph into MAGUS
# =====================================================================================================

class CSRRow:
    """Read-only dict-like view of one row of a CSR matrix (MAGUS reads its graph as graph.matrix[a][b])."""
    __slots__ = ("idx", "val")

    def __init__(self, idx, val):
        self.idx, self.val = idx, val

    def get(self, key, default=None):
        i = np.searchsorted(self.idx, key)
        return int(self.val[i]) if i < len(self.idx) and self.idx[i] == key else default

    def __getitem__(self, key):
        v = self.get(key)
        if v is None:
            raise KeyError(key)
        return v

    def __contains__(self, key):
        return self.get(key) is not None

    def items(self):
        return zip(self.idx.tolist(), self.val.tolist())

    def keys(self):
        return self.idx.tolist()

    def __iter__(self):
        return iter(self.idx.tolist())

    def __len__(self):
        return len(self.idx)


class CSRGraph:
    """MAGUS's graph.matrix (a list of dicts) replaced by a CSR matrix: same reads, ~10x less memory."""

    def __init__(self, csr):
        csr.sort_indices()
        self.csr = csr

    def __len__(self):
        return self.csr.shape[0]

    def __getitem__(self, a):
        s, e = self.csr.indptr[a], self.csr.indptr[a + 1]
        return CSRRow(self.csr.indices[s:e], self.csr.data[s:e])

    def __iter__(self):
        return (self[a] for a in range(len(self)))


def install_filter(variant, info_path):
    """Replace MAGUS's graph construction by build_evidence + filtered_graph (every merge MAGUS does)."""
    from magus.align.merge import alignment_graph
    from magus.align.merge.graph_build import graph_builder as gb
    from magus.configuration import Configs

    def build_matrix(context):
        if Configs.graphBuildRestrict:
            raise SystemExit("--graphbuildrestrict true is not supported")
        json.dump([os.path.abspath(p) for p in context.subalignmentPaths],
                  open(info_path.replace(".model.json", ".subsets.json"), "w"))
        weight, support, cover, B = build_evidence(context)
        info = {"variant": variant}
        graph = filtered_graph(context, weight, support, cover, B, variant, info).tocsr()
        graph.indices = graph.indices.astype(np.int32)
        context.graph.matrix = CSRGraph(graph)
        json.dump(info, open(info_path, "w"), indent=1, default=float)
        Configs.log("[filter_gcm] {}: kept {} of {} cross-subset edges".format(variant, info["kept_edges"],
                                                                               info["edges"]))

    def write_graph(self, path):
        """The graph file MCL reads: one 'a b weight' line per entry, rows in order, columns sorted."""
        coo = self.matrix.csr.tocoo()
        np.savetxt(path, np.column_stack([coo.row, coo.col, coo.data]), fmt="%d")

    def clustering_cost(self, clusters):
        """MAGUS's cut cost (only logged here), vectorised."""
        label = np.arange(self.matrixSize) + len(clusters)
        for i, cluster in enumerate(clusters):
            label[np.asarray(cluster, dtype=np.int64)] = i
        subset_of = np.array([s for s, _ in self.matSubPosMap])
        coo = self.matrix.csr.tocoo()
        cut = (subset_of[coo.row] != subset_of[coo.col]) & (label[coo.row] != label[coo.col])
        return int(coo.data[cut].sum() / 2)

    gb.buildMatrix = build_matrix
    alignment_graph.AlignmentGraph.writeGraphToFile = write_graph
    alignment_graph.AlignmentGraph.computeClusteringCost = clustering_cost


def magus_inner(variant, info_path, magus_args):
    """Runs inside a fresh Python process: install the filter, then run MAGUS's own main()."""
    install_filter(variant, info_path)
    from magus.main import main as magus_main
    sys.argv = ["magus"] + magus_args
    magus_main()


def run_magus(variant, magus_args, out, log_path):
    """MAGUS in a child process (MAGUS keeps global state), with the filter installed. Returns wall seconds."""
    info_path = out + ".model.json"
    start = time.time()
    with open(log_path, "w") as log:
        subprocess.run([sys.executable, os.path.abspath(__file__), "--_inner", variant, info_path] + magus_args,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    return round(time.time() - start, 1)


# =====================================================================================================
# 5. Inputs: subset order, scoring
# =====================================================================================================

def subset_order(workdir):
    """Subset alignment files in the order MAGUS numbered them in the logged run. MAGUS takes the files of a
    directory in os.listdir order, which differs between machines, and graph node ids follow that order, so the
    bank stores it in subsets.json (absolute or relative paths; only the file names are used here). Without
    subsets.json the files are sorted by subset number."""
    sub_dir = os.path.join(workdir, "inputs", "subalignments")
    names = os.listdir(sub_dir)
    order_file = os.path.join(workdir, "subsets.json")
    if os.path.exists(order_file):
        order = [os.path.basename(p) for p in json.load(open(order_file))]
        assert sorted(order) == sorted(names), "subsets.json does not match " + sub_dir
    else:
        order = sorted(names, key=lambda f: int(f.rsplit("_", 1)[1].split(".")[0]))
    return [os.path.join(sub_dir, f) for f in order]


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(seqs, path):
    with open(path, "w") as f:
        for k, v in seqs.items():
            f.write(">{}\n{}\n".format(k, v))


def fastsp(ref_path, est_path):
    """SPFN, SPFP, TC of the estimate against the reference (FastSP, upper-cased; the estimate is restricted to the
    reference's sequences with all-gap columns removed when it holds more, e.g. HomFam)."""
    ref, est = read_fasta(ref_path), read_fasta(est_path)
    if set(est) != set(ref):
        rows = [est[t] for t in ref]
        cols = [i for i in range(len(rows[0])) if any(r[i] not in "-." for r in rows)]
        est = {t: "".join(est[t][i] for i in cols) for t in ref}
    jar = os.environ.get("FASTSP_JAR", "/opt/tools/FastSP/FastSP.jar")
    with tempfile.TemporaryDirectory() as tmp:
        r, e = os.path.join(tmp, "ref.fa"), os.path.join(tmp, "est.fa")
        write_fasta({k: v.upper() for k, v in ref.items()}, r)
        write_fasta({k: v.upper() for k, v in est.items()}, e)
        out = subprocess.run(["java", "-Xmx" + os.environ.get("FASTSP_XMX", "4g"), "-jar", jar, "-r", r, "-e", e],
                             capture_output=True, text=True, check=True)
    s = {}
    for line in (out.stdout + out.stderr).splitlines():
        tok = line.split()
        if len(tok) == 2 and tok[0] in ("SPFN", "SPFP", "TC"):
            s[tok[0]] = float(tok[1])
        elif line.startswith("MaxLenNoGap"):
            fields = dict(f.strip().split("= ") for f in line.split(","))
            s["LenRef"], s["LenEst"] = int(fields["LenRef"]), int(fields["LenEst"])
    s["avgErr"] = (s["SPFN"] + s["SPFP"]) / 2
    return s


# =====================================================================================================
# 6. Command line
# =====================================================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("input", help="WORKDIR (merge only) or UNALIGNED.fasta (with --from-scratch)")
    ap.add_argument("variant")
    ap.add_argument("--out", required=True, help="merged alignment (FASTA)")
    ap.add_argument("--from-scratch", action="store_true", help="run MAGUS end to end on unaligned sequences")
    ap.add_argument("--subsets", type=int, default=25, help="number of subsets for --from-scratch (default 25)")
    ap.add_argument("--backbones", help="backbone directory (default WORKDIR/inputs/backbones)")
    ap.add_argument("--ref", help="reference alignment: score the result with FastSP and print a JSON row")
    ap.add_argument("-np", "--threads", type=int, default=4)
    ap.add_argument("--keep-work", action="store_true", help="keep MAGUS's working directory (OUT.work)")
    a = ap.parse_args()

    out = os.path.abspath(a.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    work = out + ".work"
    shutil.rmtree(work, ignore_errors=True)  # MAGUS reuses an existing graph file, so always start clean
    if a.from_scratch:
        magus_args = ["-i", os.path.abspath(a.input)] + scratch_flags(a.subsets)
    else:
        bb = os.path.abspath(a.backbones or os.path.join(a.input, "inputs", "backbones"))
        magus_args = ["-s"] + subset_order(os.path.abspath(a.input)) + ["-b", bb] + MERGE_FLAGS
    magus_args += ["-np", str(a.threads), "-d", work, "-o", out]
    wall = run_magus(a.variant, magus_args, out, out + ".log")
    if not a.keep_work:
        shutil.rmtree(work, ignore_errors=True)
    info = json.load(open(out + ".model.json"))
    row = {"variant": a.variant, "merge_wall": wall, "edges": info["edges"], "kept_edges": info["kept_edges"]}
    if a.ref:
        row.update(fastsp(a.ref, out))
    print(json.dumps(row))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--_inner":
        magus_inner(sys.argv[2], sys.argv[3], sys.argv[4:])
    else:
        main()
