"""Learned alignment-graph edge weights for GCM (MAGUS's merge).

MAGUS weights the edge between subset columns a and b by w_ab, the number of
residue pairs that the backbone alignments put in one column. Our oracle runs
show that the evidence is the main error source: with the true alignment as
the only backbone (edge weight t_ab = true homologous pairs between a and b)
MAGUS's error on 1000M2 drops from 8.2% to 4.8%. Here a regression model
predicts the fraction of truly homologous pairs y_ab = t_ab / (n_a n_b) from
cheap per-edge features of the ordinary evidence, and the merge uses the
predicted number of true pairs, round(100 * y_hat * n_a * n_b), as the weight
(the diagonal gets 100 * n_a^2, the same convention). Nothing else changes:
same subsets, same backbones, MCL, minclusters trace.

Features of an edge (a, b), a != b (same-subset edges included, as in MAGUS):
  f        w_ab / sum_k c_k(a) c_k(b)  (share of co-occurring backbone pairs aligned)
  logw     log(1 + w_ab)
  ksup     number of backbones that support the edge
  ra, rb   w_ab / max over b' in subset(b) of w_ab' (and symmetric)
  rfa, rfb same with f
  dpos     |pos_a / L_i - pos_b / L_j| (relative column positions)
  occ_a, occ_b   column occupancy n_a / s_i, n_b / s_j
  logs_i, logs_j subset sizes
  same     1 if a and b are columns of the same subset alignment

    GCMX_EDGE_DUMP=out.npz GCMX_EDGE_TRUE=true.fasta python -m gcmx.run_magus ...   # training data
    python -m gcmx.learnweights train MODEL.joblib DUMP.npz [DUMP.npz ...]
    GCMX_EDGE_MODEL=MODEL.joblib python -m gcmx.run_magus ...                        # learned weights
"""

import os
import sys

import numpy as np
import scipy.sparse as sp

from . import fasta

FEATURES = ["f", "logw", "ksup", "ra", "rb", "rfa", "rfb", "dpos", "occ_a", "occ_b", "logs_i", "logs_j", "same"]


def node_info(context):
    """Per graph node: subset index, relative position, occupancy (letters in the column), subset size."""
    graph = context.graph
    sub = np.array([s for s, _ in graph.matSubPosMap])
    pos = np.array([p for _, p in graph.matSubPosMap], dtype=float)
    lengths = np.array(graph.subalignmentLengths, dtype=float)
    occ = np.zeros(graph.matrixSize)
    size = np.zeros(len(lengths))
    for k, path in enumerate(context.subalignmentPaths):
        aln = fasta.read(path)
        size[k] = len(aln)
        counts = np.zeros(int(lengths[k]))
        for seq in aln.values():
            counts += np.frombuffer(seq.encode(), dtype=np.uint8) != ord("-")
        occ[graph.subsetMatrixIdx[k]:graph.subsetMatrixIdx[k] + int(lengths[k])] = counts
    return sub, pos / np.maximum(lengths[sub] - 1, 1), occ, size


def _group_max(key, value):
    order = np.argsort(key, kind="stable")
    k, v = key[order], value[order]
    starts = np.r_[0, np.flatnonzero(np.diff(k)) + 1]
    maxes = np.maximum.reduceat(v, starts)
    out = np.empty_like(value)
    out[order] = np.repeat(maxes, np.diff(np.r_[starts, len(k)]))
    return out


def edge_features(context, total, occupancy, support):
    """Feature matrix for the off-diagonal edges of `total` (both triangles; MAGUS keeps same-subset edges)."""
    sub, rel, occ, size = node_info(context)
    coo = total.tocoo()
    keep = coo.row != coo.col
    a, b, w = coo.row[keep], coo.col[keep], coo.data[keep].astype(float)
    C = np.vstack(occupancy)  # backbones x nodes
    denom = np.einsum("ke,ke->e", C[:, a], C[:, b])
    f = w / np.maximum(denom, 1)
    ksup = np.asarray(support[a, b]).ravel().astype(float)
    nsub = len(size)
    ra = w / _group_max(a * nsub + sub[b], w)
    rb = w / _group_max(b * nsub + sub[a], w)
    rfa = f / _group_max(a * nsub + sub[b], f)
    rfb = f / _group_max(b * nsub + sub[a], f)
    X = np.column_stack([f, np.log1p(w), ksup, ra, rb, rfa, rfb, np.abs(rel[a] - rel[b]),
                         occ[a] / size[sub[a]], occ[b] / size[sub[b]], np.log(size[sub[a]]), np.log(size[sub[b]]), (sub[a] == sub[b]).astype(float)])
    return a, b, X, occ


def true_pairs(context, true_path, n):
    from magus.helpers import sequenceutils
    from .fastgraph import _column_matrix
    from .weighting import _backbone_alignmap
    # MAGUS maps only backbone taxa to their subset-alignment rows; the true alignment has every taxon
    added = []
    for path in context.subalignmentPaths:
        for taxon, seq in sequenceutils.readFromFasta(path, removeDashes=False).items():
            if taxon not in context.backboneSubalignment:
                context.backboneSubalignment[taxon] = seq
                added.append(taxon)
    A = _column_matrix(_backbone_alignmap(context, true_path), n)
    for taxon in added:
        del context.backboneSubalignment[taxon]
    return (A.T @ A).tocsr()


def transform(context, total, occupancy, support):
    """Called by fastgraph.buildMatrix: dump training data and/or replace the edge weights."""
    from magus.configuration import Configs
    n = total.shape[0]
    a, b, X, occ = edge_features(context, total, occupancy, support)
    if os.environ.get("GCMX_EDGE_DUMP"):
        T = true_pairs(context, os.environ["GCMX_EDGE_TRUE"], n)
        t = np.asarray(T[a, b]).ravel().astype(float)
        np.savez_compressed(os.environ["GCMX_EDGE_DUMP"], X=X.astype(np.float32), t=t, nab=occ[a] * occ[b])
        Configs.log("[gcmx:learnweights] dumped {} edges".format(len(a)))
    if not os.environ.get("GCMX_EDGE_MODEL"):
        return total
    import joblib
    model = joblib.load(os.environ["GCMX_EDGE_MODEL"])
    y = np.clip(model.predict(X), 0.0, 1.0)
    w = np.rint(100 * y * occ[a] * occ[b]).astype(np.int64)
    nz = w > 0
    diag = np.arange(n)
    rows = np.r_[a[nz], diag]
    cols = np.r_[b[nz], diag]
    data = np.r_[w[nz], np.rint(100 * occ * occ).astype(np.int64)]
    Configs.log("[gcmx:learnweights] reweighted {} edges ({} kept)".format(len(a), int(nz.sum())))
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n))


def train(model_path, dumps, max_edges=1_500_000, seed=0):
    import joblib
    from sklearn.ensemble import HistGradientBoostingRegressor
    rng = np.random.default_rng(seed)
    Xs, ys, ws = [], [], []
    per = max_edges // len(dumps)
    for path in dumps:
        d = np.load(path)
        X, t, nab = d["X"], d["t"], d["nab"]
        idx = rng.choice(len(t), size=min(per, len(t)), replace=False)
        Xs.append(X[idx]); ys.append(np.clip(t[idx] / np.maximum(nab[idx], 1), 0, 1)); ws.append(nab[idx])
    X, y, w = np.vstack(Xs), np.concatenate(ys), np.concatenate(ws)
    model = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.1, max_leaf_nodes=63, early_stopping=True,
                                          validation_fraction=0.1, n_iter_no_change=10, random_state=seed)
    model.fit(X, y, sample_weight=w)
    joblib.dump(model, model_path)
    pred = np.clip(model.predict(X), 0, 1)
    base = np.clip(X[:, 0], 0, 1)
    err = lambda p: float(np.average((p - y) ** 2, weights=w))
    print("train edges {}  weighted MSE: model {:.4f}  vs raw f {:.4f}".format(len(y), err(pred), err(base)))


if __name__ == "__main__":
    if sys.argv[1] == "train":
        train(sys.argv[2], sys.argv[3:])
