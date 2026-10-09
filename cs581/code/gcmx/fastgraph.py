"""Vectorized, memory-compact alignment graph (drop-in for MAGUS's buildMatrix).

MAGUS adds, for every backbone column l and every pair of graph nodes (a, b)
with letters in l, count_l(a) * count_l(b) to edge (a, b), in a Python double
loop, and stores the result as a list of Python dicts. With A_b the sparse
(backbone columns x nodes) count matrix of backbone b, the graph is exactly
sum_b A_b^T A_b, which scipy computes in C. Self-loops and same-subalignment
edges are kept, as in MAGUS (--graphbuildrestrict is honoured).

The result is stored as a CSR matrix behind a read-only, dict-like interface
(graph.matrix[a].items() / .get(b, 0) / [b] / `in`), which is all MAGUS's
clustering, cleanup, trace and optimizer code uses. With HMM-extended
backbones and soft constraints the graph reaches ~50M entries, which as
Python dicts needs >7 GB; as CSR it needs ~0.6 GB. Graph-file writing and the
clustering cost are vectorized as well. Weights stay integers, so MCL sees the
same graph.
"""

import json
import os

import numpy as np
import scipy.sparse as sp

from magus.align.merge import alignment_graph
from magus.align.merge.graph_build import graph_builder as gb
from magus.configuration import Configs
from magus.tasks import task

from . import fasta
from .weighting import _backbone_alignmap


class Row:
    """Read-only dict-like view of one CSR row (sorted column indices)."""
    __slots__ = ("idx", "val")

    def __init__(self, idx, val):
        self.idx, self.val = idx, val

    def get(self, key, default=None):
        i = np.searchsorted(self.idx, key)
        if i < len(self.idx) and self.idx[i] == key:
            return int(self.val[i])
        return default

    def __getitem__(self, key):
        value = self.get(key)
        if value is None:
            raise KeyError(key)
        return value

    def __contains__(self, key):
        return self.get(key) is not None

    def items(self):
        return zip(self.idx.tolist(), self.val.tolist())

    def keys(self):
        return self.idx.tolist()

    def values(self):
        return self.val.tolist()

    def __iter__(self):
        return iter(self.idx.tolist())

    def __len__(self):
        return len(self.idx)


class CSRGraph:
    def __init__(self, csr):
        csr.sort_indices()
        self.csr = csr

    def __len__(self):
        return self.csr.shape[0]

    def __getitem__(self, a):
        s, e = self.csr.indptr[a], self.csr.indptr[a + 1]
        return Row(self.csr.indices[s:e], self.csr.data[s:e])

    def __iter__(self):
        return (self[a] for a in range(len(self)))


def install():
    gb.buildMatrix = buildMatrix
    original_write = alignment_graph.AlignmentGraph.writeGraphToFile
    original_cost = alignment_graph.AlignmentGraph.computeClusteringCost

    def writeGraphToFile(self, filePath):
        if not isinstance(self.matrix, CSRGraph):
            return original_write(self, filePath)
        coo = self.matrix.csr.tocoo()
        np.savetxt(filePath, np.column_stack([coo.row, coo.col, coo.data]), fmt="%d")
        Configs.log("Wrote matrix to {}".format(filePath))

    def computeClusteringCost(self, clusters):
        if not isinstance(self.matrix, CSRGraph):
            return original_cost(self, clusters)
        label = np.arange(self.matrixSize) + len(clusters)  # unclustered nodes: own cluster
        for n, cluster in enumerate(clusters):
            label[np.asarray(cluster, dtype=np.int64)] = n
        sub = np.array([s for s, _ in self.matSubPosMap])
        coo = self.matrix.csr.tocoo()
        cut = (sub[coo.row] != sub[coo.col]) & (label[coo.row] != label[coo.col])
        return int(coo.data[cut].sum() / 2)

    alignment_graph.AlignmentGraph.writeGraphToFile = writeGraphToFile
    alignment_graph.AlignmentGraph.computeClusteringCost = computeClusteringCost


def _column_matrix(alignmap, n):
    rows, cols, vals = [], [], []
    for l, column in enumerate(alignmap):
        for node, count in column.items():
            rows.append(l)
            cols.append(node)
            vals.append(count)
    return sp.csr_matrix((np.array(vals, dtype=np.int64), (rows, cols)), shape=(len(alignmap), n))


PP_VALUE = {"*": 0.975, **{str(d): (0.025 if d == 0 else d / 10.0) for d in range(10)}}


def _weighted_column_matrix(context, alignedFile, pp_path, n):
    """Like MAGUS's backboneToAlignMap, but each residue counts with its HMMER posterior
    (residues without a posterior, i.e. the backbone's own sequences, count 1)."""
    graph = context.graph
    backbone = fasta.read(alignedFile)
    pps = fasta.read(pp_path)
    rows, cols, vals = [], [], []
    for taxon, backboneseq in backbone.items():
        sub = context.taxonSubalignmentMap[taxon]
        subsetseq = context.backboneSubalignment[taxon].seq
        unaligned = context.unalignedSequences[taxon].seq
        posarray = [n_ for n_, ch in enumerate(subsetseq) if ch not in "-."]
        weights = [PP_VALUE.get(ch, 1.0) for ch in pps[taxon]] if taxon in pps else None
        i = col = 0
        for c in backboneseq:
            if i == len(unaligned):
                break
            if c == unaligned[i]:
                rows.append(col)
                cols.append(graph.subsetMatrixIdx[sub] + posarray[i])
                vals.append(weights[i] if weights else 1.0)
            if c.upper() == unaligned[i]:
                i += 1
            if c == c.upper() and c != ".":
                col += 1
    length = max(rows) + 1 if rows else 1
    return sp.csr_matrix((np.array(vals), (rows, cols)), shape=(length, n))


def buildMatrix(context):
    graph = context.graph
    n = graph.matrixSize
    files = [t.outputFile for t in task.asCompleted(context.backboneTasks)]
    for path in context.backbonePaths:
        if path not in files:
            files.append(path)

    # optional per-backbone weights: GCMX_BACKBONE_WEIGHTS = JSON file {basename: weight};
    # weights are scaled by 100 and rounded so that graph weights stay integers
    weights = {}
    if os.environ.get("GCMX_BACKBONE_WEIGHTS"):
        weights = json.load(open(os.environ["GCMX_BACKBONE_WEIGHTS"]))
    total = sp.csr_matrix((n, n), dtype=np.int64)
    learn = bool(os.environ.get("GCMX_EDGE_DUMP") or os.environ.get("GCMX_EDGE_MODEL"))
    occupancy, support = [], sp.csr_matrix((n, n), dtype=np.int64)  # for gcmx.learnweights
    for alignedFile in files:
        Configs.log("[gcmx:fastgraph] Feeding backbone {} to the graph..".format(alignedFile))
        pp_path = os.path.join(os.environ.get("GCMX_PP_DIR", ""), os.path.basename(alignedFile))
        if os.environ.get("GCMX_PP_DIR") and os.path.exists(pp_path):
            A = _weighted_column_matrix(context, alignedFile, pp_path, n)
            product = (A.T @ A).tocsr()
            product = sp.csr_matrix((np.rint(100 * product.data).astype(np.int64), product.indices, product.indptr),
                                    shape=product.shape)  # keep integer weights (x100) for MCL
        else:
            A = _column_matrix(_backbone_alignmap(context, alignedFile), n)
            product = (A.T @ A).tocsr()
            if learn:
                occupancy.append(np.asarray(A.sum(axis=0)).ravel().astype(float))
                support = support + (product > 0).astype(np.int64)
            if os.environ.get("GCMX_PP_DIR"):
                product = product * 100  # same scale as posterior-weighted backbones
        if weights:
            product = product * int(round(100 * weights.get(os.path.basename(alignedFile), 1.0)))
        total = total + product

    if learn:
        from . import learnweights
        total = learnweights.transform(context, total.tocsr(), occupancy, support.tocsr())

    if Configs.graphBuildRestrict:
        coo = total.tocoo()
        sub = np.array([s for s, _ in graph.matSubPosMap])
        keep = (sub[coo.row] != sub[coo.col]) | (coo.row == coo.col)
        total = sp.csr_matrix((coo.data[keep], (coo.row[keep], coo.col[keep])), shape=(n, n))
    total.indices = total.indices.astype(np.int32)
    graph.matrix = CSRGraph(total)
    Configs.log("[gcmx:fastgraph] Built graph with {} weighted entries from {} backbones".format(total.nnz, len(files)))
