"""GCM merge where each backbone file's residue-pair counts are multiplied by an integer weight.

    GG_WEIGHTS=weights.json python run_wmerge.py <gcmx.run_magus args>

weights.json maps backbone file basename -> integer weight (default 1). The graph is otherwise built
exactly as MAGUS's graph_builder.addAlignmentFileToGraph (same loop as protcons/code/run_edgesup.py), so
weight 1 everywhere is MAGUS itself and weight k equals k copies of the file.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "code")))

from magus.align.merge.graph_build import graph_builder as gb  # noqa: E402
from magus.configuration import Configs  # noqa: E402
from magus.tasks import task  # noqa: E402

from gcmx import weighting  # noqa: E402

W = json.load(open(os.environ["GG_WEIGHTS"]))
K = int(os.environ.get("GG_ESK", "1"))  # also delete edges that fewer than K of the weight-1 files support


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
        w = int(W.get(os.path.basename(alignedFile), 1))
        alignmap = weighting._backbone_alignmap(context, alignedFile)
        seen = set()
        for column in alignmap:
            items = list(column.items())
            for a, avalue in items:
                row = graph.matrix[a]
                for b, bvalue in items:
                    if Configs.graphBuildRestrict:
                        asub, apos = graph.matSubPosMap[a]
                        bsub, bpos = graph.matSubPosMap[b]
                        if asub == bsub and apos != bpos:
                            continue
                    row[b] = row.get(b, 0) + w * avalue * bvalue
                    if w == 1:
                        seen.add((a, b))
        for a, b in seen:
            nbb[a][b] = nbb[a].get(b, 0) + 1
    if K > 1:
        removed = 0
        for a in range(graph.matrixSize):
            row = graph.matrix[a]
            for b in list(row):
                if a != b and nbb[a].get(b, 0) < K:
                    del row[b]
                    removed += 1
        Configs.log("[wmerge] K={}: removed {} directed edges".format(K, removed))
    Configs.log("[wmerge] weights {}".format(sorted(set(W.values()))))


gb.buildMatrix = buildMatrix
from gcmx import run_magus  # noqa: E402

run_magus.main()
