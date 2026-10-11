"""GCM merge with a global edge-support threshold (post-hoc baseline, see REPORT.md).

    python run_edgesup.py K <gcmx.run_magus args>

The alignment graph is built exactly as MAGUS does (edge weight = sum over backbones of aligned residue
pairs between the two subset columns), then every edge between two different subset columns that fewer
than K backbones contribute to is deleted. K = 1 is MAGUS itself.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "code")))

from magus.align.merge.graph_build import graph_builder as gb  # noqa: E402
from magus.configuration import Configs  # noqa: E402
from magus.tasks import task  # noqa: E402

from gcmx import weighting  # noqa: E402

K = int(sys.argv.pop(1))


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
                    if Configs.graphBuildRestrict:
                        asub, apos = graph.matSubPosMap[a]
                        bsub, bpos = graph.matSubPosMap[b]
                        if asub == bsub and apos != bpos:
                            continue
                    row[b] = row.get(b, 0) + avalue * bvalue
                    seen.add((a, b))
        for a, b in seen:
            nbb[a][b] = nbb[a].get(b, 0) + 1
    removed = kept = 0
    for a in range(graph.matrixSize):
        row = graph.matrix[a]
        for b in list(row):
            if a != b and nbb[a].get(b, 0) < K:
                del row[b]
                removed += 1
            else:
                kept += 1
    Configs.log("[edgesup] K={}: removed {} directed edges, kept {}".format(K, removed, kept))


gb.buildMatrix = buildMatrix
from gcmx import run_magus  # noqa: E402

run_magus.main()
