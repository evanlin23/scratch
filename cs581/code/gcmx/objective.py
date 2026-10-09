"""Evaluate MWT-AM objective of several traces under one alignment graph.

    python -m gcmx.objective GRAPH.txt TRACE.txt [TRACE.txt ...]

GRAPH.txt is a MAGUS graph file ("a b w" lines, nodes = constraint columns);
each TRACE.txt is a MAGUS trace (one cluster of nodes per line). All runs
must use the same subalignment files (same node numbering). Reports, per
trace, the trace (MWT) score = total weight of edges inside clusters, and the
cut = total weight between clusters (MAGUS's "total cost"), excluding edges
between columns of the same constraint alignment is not possible here, so
every non-self-loop edge counts -- compare traces only relative to each other.
"""

import sys


def load_graph(path):
    edges = []
    with open(path) as f:
        for line in f:
            a, b, w = line.split()
            a, b = int(a), int(b)
            if a < b:
                edges.append((a, b, float(w)))
    return edges


def load_trace(path):
    cluster_of = {}
    with open(path) as f:
        for i, line in enumerate(f):
            for node in line.split():
                cluster_of[int(node)] = i
    return cluster_of


def main():
    edges = load_graph(sys.argv[1])
    total = sum(w for _, _, w in edges)
    for path in sys.argv[2:]:
        cluster_of = load_trace(path)
        inside = sum(w for a, b, w in edges
                     if a in cluster_of and cluster_of[a] == cluster_of.get(b, -1))
        print("{:>14.0f} inside ({:.2%} of edge weight)  {}".format(inside, inside / total, path))


if __name__ == "__main__":
    main()
