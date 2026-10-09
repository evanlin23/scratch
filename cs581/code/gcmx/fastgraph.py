"""Vectorized alignment-graph construction (drop-in for MAGUS's buildMatrix).

MAGUS adds, for every backbone column l and every pair of graph nodes (a, b)
with letters in l, count_l(a) * count_l(b) to edge (a, b), in a Python double
loop. With A_b the sparse (backbone columns x nodes) count matrix of backbone b,
that is exactly the sum over backbones of A_b^T A_b, which scipy computes in
C. Self-loops and same-subalignment edges are kept, as in MAGUS (unless
--graphbuildrestrict is set, which this module also honours). Weights stay
integers, so the graph file and everything downstream are unchanged.
"""

import numpy as np
import scipy.sparse as sp

from magus.align.merge.graph_build import graph_builder as gb
from magus.configuration import Configs
from magus.tasks import task

from .weighting import _backbone_alignmap


def install():
    gb.buildMatrix = buildMatrix


def _column_matrix(alignmap, n):
    rows, cols, vals = [], [], []
    for l, column in enumerate(alignmap):
        for node, count in column.items():
            rows.append(l)
            cols.append(node)
            vals.append(count)
    return sp.csr_matrix((np.array(vals, dtype=np.int64), (rows, cols)), shape=(len(alignmap), n))


def buildMatrix(context):
    graph = context.graph
    n = graph.matrixSize
    files = [t.outputFile for t in task.asCompleted(context.backboneTasks)]
    for path in context.backbonePaths:
        if path not in files:
            files.append(path)

    total = sp.csr_matrix((n, n), dtype=np.int64)
    for alignedFile in files:
        Configs.log("[gcmx:fastgraph] Feeding backbone {} to the graph..".format(alignedFile))
        A = _column_matrix(_backbone_alignmap(context, alignedFile), n)
        total = total + (A.T @ A).tocsr()

    coo = total.tocoo()
    if Configs.graphBuildRestrict:
        sub = np.array([graph.matSubPosMap[i][0] for i in range(n)])
        keep = (sub[coo.row] != sub[coo.col]) | (coo.row == coo.col)
        coo = sp.coo_matrix((coo.data[keep], (coo.row[keep], coo.col[keep])), shape=(n, n))
    matrix = [dict() for _ in range(n)]
    for a, b, w in zip(coo.row.tolist(), coo.col.tolist(), coo.data.tolist()):
        matrix[a][b] = w
    graph.matrix = matrix
    Configs.log("[gcmx:fastgraph] Built graph with {} weighted entries from {} backbones".format(coo.nnz, len(files)))
