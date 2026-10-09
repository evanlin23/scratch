"""Alternative edge-weighting schemes for the GCM alignment graph.

Notation: nodes are constraint-alignment columns. For backbone b, let
  c_b(u, v) = number of letter pairs (x in column u, y in column v) that the
              backbone puts in the same column   (MAGUS default edge weight)
  o_b(u)    = number of letters of column u that appear in backbone b
              (i.e. how many backbone sequences have a letter in column u)
so o_b(u) * o_b(v) is the number of *opportunities* backbone b had to
support the pair (u, v), and c_b(u, v) / (o_b(u) o_b(v)) is in [0, 1].

Schemes:
  count       sum_b c_b(u,v)                                   (MAGUS default)
  frac        sum_b c_b(u,v) / sum_b o_b(u) o_b(v)             (pooled support fraction)
  fracshrink  sum_b c_b(u,v) / (sum_b o_b(u) o_b(v) + alpha)   (shrunk toward 0)
  vote        sum_b c_b(u,v) / (o_b(u) o_b(v))                 (each backbone = 1 vote)

The default scheme lets densely populated columns dominate: a pair of
gap-free columns can collect up to (L/k)^2 support per backbone, while a pair
of gappy columns with one sampled letter each collects at most 1. The other
schemes measure *how consistently* the backbones align two columns instead.
"""

import numpy as np

from magus.align.merge.graph_build import graph_builder as gb
from magus.configuration import Configs
from magus.helpers import sequenceutils
from magus.tasks import task

SCHEMES = ("count", "frac", "fracshrink", "vote")
_state = {"scheme": "count", "alpha": 1.0}


def install(scheme="count", alpha=1.0):
    if scheme not in SCHEMES:
        raise ValueError("unknown weighting scheme {}; choose from {}".format(scheme, SCHEMES))
    _state["scheme"] = scheme
    _state["alpha"] = alpha
    if scheme != "count":
        gb.buildMatrix = buildMatrix


def _backbone_alignmap(context, alignedFile):
    """Same front half as magus graph_builder.addAlignmentFileToGraph."""
    backboneAlign = sequenceutils.readFromFasta(alignedFile)
    alignmentLength = len(next(iter(backboneAlign.values())).seq)
    if alignedFile in context.backboneExtend:
        extensionTasks = gb.requestHmmExtensionTasks(context, backboneAlign, alignedFile)
        task.submitTasks(extensionTasks)
        for extensionTask in task.asCompleted(extensionTasks):
            backboneAlign.update(sequenceutils.readFromStockholm(extensionTask.outputFile, includeInsertions=True))
    return gb.backboneToAlignMap(context, backboneAlign, alignmentLength)


def _backbone_counts(graph, alignmap):
    """Return (c_b as dict-of-dicts, o_b as array) for one backbone."""
    occ = np.zeros(graph.matrixSize)
    counts = {}
    for column in alignmap:
        items = list(column.items())
        for a, avalue in items:
            occ[a] += avalue
            row = counts.setdefault(a, {})
            for b, bvalue in items:
                if Configs.graphBuildRestrict:
                    asub, apos = graph.matSubPosMap[a]
                    bsub, bpos = graph.matSubPosMap[b]
                    if asub == bsub and apos != bpos:
                        continue
                row[b] = row.get(b, 0) + avalue * bvalue
    return counts, occ


def buildMatrix(context):
    scheme, alpha = _state["scheme"], _state["alpha"]
    graph = context.graph

    files = [backboneTask.outputFile for backboneTask in task.asCompleted(context.backboneTasks)]
    for backboneFile in context.backbonePaths:
        if backboneFile not in files:
            files.append(backboneFile)

    weights = [dict() for _ in range(graph.matrixSize)]
    occupancies = []
    for alignedFile in files:
        Configs.log("[gcmx:{}] Feeding backbone {} to the graph..".format(scheme, alignedFile))
        counts, occ = _backbone_counts(graph, _backbone_alignmap(context, alignedFile))
        occupancies.append(occ)
        for a, row in counts.items():
            target = weights[a]
            if scheme == "vote":
                for b, value in row.items():
                    target[b] = target.get(b, 0.0) + value / (occ[a] * occ[b])
            else:
                for b, value in row.items():
                    target[b] = target.get(b, 0.0) + value

    if scheme in ("frac", "fracshrink"):
        occ = np.array(occupancies)  # backbones x nodes
        pad = alpha if scheme == "fracshrink" else 0.0
        for a in range(graph.matrixSize):
            if not weights[a]:
                continue
            nbrs = np.fromiter(weights[a].keys(), dtype=np.int64)
            vals = np.fromiter(weights[a].values(), dtype=float)
            opportunities = (occ[:, a][:, None] * occ[:, nbrs]).sum(axis=0)
            weights[a] = dict(zip(nbrs.tolist(), (vals / (opportunities + pad)).tolist()))

    graph.matrix = weights
    Configs.log("[gcmx:{}] Built weighted matrix from {} backbones".format(scheme, len(files)))
