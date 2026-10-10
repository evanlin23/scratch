"""MAGUS merge (paper flags, MCL -I 4) that also dumps the alignment graph with per-edge backbone support.

    python build_graph.py OUT.npz <gcmx.run_magus args>

The graph is built exactly as MAGUS does (edge weight = residue-pair count summed over backbones; a copy of
cs581/protcons/code/run_edgesup.py with K = 1). OUT.npz holds, for every directed edge (a, b) of the graph:
a, b, w (weight), nbb (number of backbones that contribute to the edge), plus node -> subset index.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "code")))

from magus.align.merge.graph_build import graph_builder as gb  # noqa: E402
from magus.configuration import Configs  # noqa: E402
from magus.tasks import task  # noqa: E402

from gcmx import weighting  # noqa: E402

OUT = sys.argv.pop(1)


def buildMatrix(context):
    graph = context.graph
    files = []
    for backboneTask in task.asCompleted(context.backboneTasks):
        files.append(backboneTask.outputFile)
    for backboneFile in context.backbonePaths:
        if backboneFile not in files:
            files.append(backboneFile)
    nbb = [dict() for _ in range(graph.matrixSize)]
    for alignedFile in files:
        alignmap = weighting._backbone_alignmap(context, alignedFile)
        seen = set()
        for column in alignmap:
            items = list(column.items())
            for a, avalue in items:
                row = graph.matrix[a]
                for b, bvalue in items:
                    row[b] = row.get(b, 0) + avalue * bvalue
                    seen.add((a, b))
        for a, b in seen:
            nbb[a][b] = nbb[a].get(b, 0) + 1
    A, B, W, N = [], [], [], []
    for a in range(graph.matrixSize):
        for b, w in graph.matrix[a].items():
            A.append(a); B.append(b); W.append(w); N.append(nbb[a].get(b, 0))
    sub = np.array([s for s, _ in graph.matSubPosMap], dtype=np.int32)
    pos = np.array([p for _, p in graph.matSubPosMap], dtype=np.int32)
    np.savez_compressed(OUT, a=np.array(A, np.int32), b=np.array(B, np.int32), w=np.array(W, np.int32),
                        nbb=np.array(N, np.int16), sub=sub, pos=pos, nfiles=len(files))
    Configs.log("[build_graph] dumped {} directed edges from {} backbones".format(len(A), len(files)))


gb.buildMatrix = buildMatrix
from gcmx import run_magus  # noqa: E402

run_magus.main()
