"""Progressive exact pairwise merging on the GCM alignment graph ("progdp").

Replaces MCL clustering + A* ordering. Starting from the k constraint
alignments, repeatedly pick two (merged) alignments A and B and merge them
with the *optimal* non-crossing column matching, i.e. the exact Maximum
Weight Trace for two alignments:

    s(i, j)  = sum of graph weights between the nodes in column a_i and b_j
    D[i][j]  = max(D[i-1][j], D[i][j-1], D[i-1][j-1] + s(i, j))   (s > 0 only)

with zero gap penalty (as in WITCH-NG's two-alignment MWT). Each row is
vectorised: D[i] = cummax(max(D[i-1][1:], D[i-1][:-1] + s_i)).

Optional refinement ("--gcmx-refine R"): R rounds of leave-one-out
re-merging. Each constraint alignment is removed from the merged alignment
and merged back with the exact two-alignment DP. Its old placement is a
feasible solution of that DP, so the total MWT score never decreases
(coordinate ascent on MWT-AM with exact steps).

Merge order ("--gcmx-order"):
  upgma  merge the pair with the highest average inter-alignment weight
         T(A,B) / (|A| |B|) (|.| = number of constraint alignments), then
         update T by summation (UPGMA on similarities)
  grow   start from the most similar pair, then repeatedly add the
         constraint alignment with the largest total weight to the growing
         alignment
"""

import numpy as np
import scipy.sparse as sp

from magus.align.merge import merger
from magus.align.merge.graph_trace import tracer
from magus.configuration import Configs

_state = {"order": "upgma", "refine": 0}


def install(order="upgma", refine=0):
    _state["order"] = order
    _state["refine"] = refine
    original = tracer.findTrace

    def findTrace(graph):
        if Configs.graphTraceMethod != "progdp":
            return original(graph)
        progressiveMerge(graph, _state["order"], _state["refine"])
        graph.writeClustersToFile(graph.tracePath)
        Configs.log("Found a trace with {} clusters and a total cost of {}".format(
            len(graph.clusters), graph.computeClusteringCost(graph.clusters)))

    tracer.findTrace = findTrace
    merger.findTrace = findTrace


def interSubalignmentWeights(graph):
    csr = getattr(graph.matrix, "csr", None)  # fastgraph's compact graph
    if csr is not None:
        coo = csr.tocoo()
        sub = np.array([s for s, _ in graph.matSubPosMap])
        keep = sub[coo.row] != sub[coo.col]
        return sp.csr_matrix((coo.data[keep].astype(float), (coo.row[keep], coo.col[keep])), shape=csr.shape)
    rows, cols, vals = [], [], []
    for a in range(graph.matrixSize):
        asub = graph.matSubPosMap[a][0]
        for b, value in graph.matrix[a].items():
            if graph.matSubPosMap[b][0] != asub:
                rows.append(a)
                cols.append(b)
                vals.append(value)
    n = graph.matrixSize
    return sp.csr_matrix((np.array(vals, dtype=float), (rows, cols)), shape=(n, n))


def membership(columns, n):
    rows = [i for i, column in enumerate(columns) for _ in column]
    nodes = [node for column in columns for node in column]
    return sp.csr_matrix((np.ones(len(nodes)), (rows, nodes)), shape=(len(columns), n))


def mergePair(A, B, W):
    """Exact two-alignment MWT merge; A, B are ordered lists of node lists."""
    n = W.shape[0]
    S = (membership(A, n) @ (W @ membership(B, n).T)).toarray()
    p, q = S.shape
    D = np.zeros((p + 1, q + 1))
    for i in range(1, p + 1):
        s = S[i - 1]
        diag = np.where(s > 0, D[i - 1, :-1] + s, -np.inf)
        D[i, 1:] = np.maximum.accumulate(np.maximum(D[i - 1, 1:], diag))

    merged = []
    i, j = p, q
    while i > 0 and j > 0:
        if D[i, j] == D[i - 1, j]:
            merged.append(A[i - 1])
            i -= 1
        elif D[i, j] == D[i, j - 1]:
            merged.append(B[j - 1])
            j -= 1
        else:
            merged.append(A[i - 1] + B[j - 1])
            i -= 1
            j -= 1
    merged.extend(reversed(A[:i]))
    merged.extend(reversed(B[:j]))
    merged.reverse()
    return merged, D[p, q]


def refine(columns, graph, W, rounds):
    k = len(graph.subalignmentLengths)
    sub = {}
    for i in range(k):
        for j in range(graph.subalignmentLengths[i]):
            sub[graph.subsetMatrixIdx[i] + j] = i
    for r in range(rounds):
        changed = 0
        for i in range(k):
            rest = [[n for n in column if sub[n] != i] for column in columns]
            rest = [column for column in rest if column]
            own = [[graph.subsetMatrixIdx[i] + j] for j in range(graph.subalignmentLengths[i])]
            before = set(tuple(sorted(column)) for column in columns)
            columns, _ = mergePair(rest, own, W)
            changed += len(set(tuple(sorted(column)) for column in columns) - before)
        Configs.log("[gcmx:progdp] refinement round {}: {} columns changed".format(r + 1, changed))
        if changed == 0:
            break
    return columns


def progressiveMerge(graph, order, rounds=0):
    k = len(graph.subalignmentLengths)
    W = interSubalignmentWeights(graph)
    Configs.log("[gcmx:progdp] {} constraint alignments, {} weighted inter-alignment edges, order={}".format(
        k, W.nnz, order))

    starts = graph.subsetMatrixIdx
    groups = {i: [[starts[i] + j] for j in range(graph.subalignmentLengths[i])] for i in range(k)}
    sizes = {i: 1 for i in range(k)}

    # T[i][j] = total edge weight between constraint alignments i and j
    sub = np.zeros(graph.matrixSize, dtype=int)
    for i in range(k):
        sub[starts[i]:starts[i] + graph.subalignmentLengths[i]] = i
    coo = W.tocoo()
    T = np.zeros((k, k))
    np.add.at(T, (sub[coo.row], sub[coo.col]), coo.data)
    T = {i: {j: T[i, j] for j in range(k) if j != i} for i in range(k)}

    def similarity(a, b):
        return T[a][b] / (sizes[a] * sizes[b])

    total = 0.0
    if order == "grow":
        a, b = max(((a, b) for a in T for b in T[a] if a < b), key=lambda ab: T[ab[0]][ab[1]])
        current, score = mergePair(groups.pop(a), groups.pop(b), W)
        total += score
        support = {c: T[a][c] + T[b][c] for c in groups}
        while groups:
            c = max(support, key=support.get)
            current, score = mergePair(current, groups.pop(c), W)
            total += score
            support.pop(c)
            for d in support:
                support[d] += T[c][d]
        final = current
    else:
        while len(groups) > 1:
            a, b = max(((a, b) for a in groups for b in groups if a < b), key=lambda ab: similarity(*ab))
            merged, score = mergePair(groups.pop(a), groups.pop(b), W)
            total += score
            new = max(T) + 1
            T[new] = {}
            for c in groups:
                T[new][c] = T[c][new] = T[a][c] + T[b][c]
            for x in (a, b):
                for c in T.pop(x):
                    T.get(c, {}).pop(x, None)
            groups[new] = merged
            sizes[new] = sizes.pop(a) + sizes.pop(b)
        final = next(iter(groups.values()))

    if rounds:
        final = refine(final, graph, W, rounds)
    graph.clusters = [sorted(column) for column in final]
    Configs.log("[gcmx:progdp] merged into {} columns; sum of pairwise-merge MWT scores {}".format(
        len(final), total))
