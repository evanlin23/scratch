# AI-assisted (Claude), exploration code for CS581 project
"""Where do the false homology pairs of a GCM-merged alignment come from?

    python3 fpsource.py KEY VARIANT [...]  -> data/fpsource.jsonl

Every pair of residues in one final column is either
  within-subset   both residues from the same subset: aligned by MAGUS's subset alignment (merge can't change it)
  direct edge     cross-subset, and the GCM graph has an edge between the two nodes, split by its posterior
                  (beta-binomial vote model, post_bb) into bins
  transitive      cross-subset, no edge between the two nodes: put together by clustering / ordering
Edge classes are further split by whether the two nodes' majority true columns agree ([edge right]: the
false pairs come from subset columns that already mix true columns) or not ([edge wrong]).
For each class: true pairs, false pairs (vs the true alignment), as % of all estimated pairs. Also, for the
cherries' share: residue pairs between true-tree cherry taxa are reported separately elsewhere (measures2).
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from measures import read_fasta, W  # noqa: E402

BINS = [0, 0.5, 0.9, 0.99, 1.0001]


def run(key, variant):
    rep = os.path.join(W, key)
    vd = os.path.join(rep, "vote", variant)
    aln = read_fasta(os.path.join(vd, "out.fasta"))
    true = read_fasta(os.path.join(rep, "true.fasta"))
    E = np.load(os.path.join(vd, "edges.npz"))
    order = json.load(open(os.path.join(vd, "subsets.json")))
    # per node: subset, final column, histogram of true columns
    node_sub, node_final, node_hist = [], [], []
    for si, path in enumerate(order):
        sa = read_fasta(path)
        L = len(next(iter(sa.values())))
        fin = np.full(L, -1, dtype=np.int64)
        hist = [defaultdict(int) for _ in range(L)]
        for t, s in sa.items():
            sc = [i for i, ch in enumerate(s) if ch not in "-."]
            fc = [i for i, ch in enumerate(aln[t]) if ch not in "-."]
            tc = [i for i, ch in enumerate(true[t]) if ch not in "-."]
            for x, y, z in zip(sc, fc, tc):
                fin[x] = y
                hist[x][z] += 1
        node_sub += [si] * L
        node_final += list(fin)
        node_hist += hist
    node_sub, node_final = np.array(node_sub), np.array(node_final)
    edge = {}
    for a, b, p in zip(E["a"], E["b"], E["post_bb"]):
        edge[(min(a, b), max(a, b))] = p
    cols = defaultdict(list)
    for v, c in enumerate(node_final):
        if c >= 0:
            cols[c].append(v)
    acc = defaultdict(lambda: [0, 0])  # class -> [true pairs, false pairs]
    for c, nodes in cols.items():
        for v in nodes:  # within one node: same subset column (subset alignment)
            h = node_hist[v]
            tot = sum(h.values())
            tp = sum(x * (x - 1) // 2 for x in h.values())
            acc["within-subset"][0] += tp
            acc["within-subset"][1] += tot * (tot - 1) // 2 - tp
        for i in range(len(nodes)):
            hi = node_hist[nodes[i]]
            ni = sum(hi.values())
            for j in range(i + 1, len(nodes)):
                hj = node_hist[nodes[j]]
                nj = sum(hj.values())
                tp = sum(x * hj.get(k, 0) for k, x in hi.items())
                u, w = nodes[i], nodes[j]
                p = edge.get((min(u, w), max(u, w)))
                if p is None:
                    cl = "transitive (no edge)"
                else:
                    k = np.searchsorted(BINS, p, side="right") - 1
                    cl = "edge post {}-{}".format(BINS[k], min(BINS[k + 1], 1))
                    # node-level: is the edge right (same majority true column) or wrong?
                    mi, mj = max(hi, key=hi.get), max(hj, key=hj.get)
                    cl += " [edge right, impure nodes]" if mi == mj else " [edge wrong]"
                acc[cl][0] += tp
                acc[cl][1] += ni * nj - tp
    tot = sum(a + b for a, b in acc.values())
    fp = sum(b for a, b in acc.values())
    row = {"key": key, "variant": variant, "est_pairs": tot, "SPFP": fp / tot,
           "classes": {k: {"true_pct_of_est": 100 * a / tot, "false_pct_of_est": 100 * b / tot,
                           "share_of_false": 100 * b / fp, "precision": a / max(a + b, 1)} for k, (a, b) in acc.items()}}
    return row


if __name__ == "__main__":
    args = sys.argv[1:]
    with open(os.path.join(HERE, "..", "data", "fpsource.jsonl"), "a") as f:
        for i in range(0, len(args), 2):
            r = run(args[i], args[i + 1])
            f.write(json.dumps(r) + "\n")
            print(r["key"], r["variant"], "SPFP %.2f" % (100 * r["SPFP"]))
            for k, v in sorted(r["classes"].items()):
                print("   {:24s} share of FP {:5.1f}%  precision {:.3f}  (% of est pairs: true {:.2f} false {:.2f})".format(
                    k, v["share_of_false"], v["precision"], v["true_pct_of_est"], v["false_pct_of_est"]))
