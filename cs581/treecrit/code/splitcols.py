# AI-assisted (Claude), exploration code for CS581 project
"""Posterior-guided column splitting of a GCM-merged alignment (Part B candidate).

    python3 splitcols.py VOTEDIR TAU [--out OUT.fasta] [--post post_bb|post_binom|post]

VOTEDIR is a vote.py output directory (out.fasta, edges.npz, subsets.json). Every final column is a set of
nodes (one column of one subset alignment each). Inside a final column, link two nodes when the cross-subset
GCM edge between them has posterior > TAU (beta-binomial vote model by default). Each connected component
becomes its own column (components of a column are written next to each other, so every row keeps its
residue order). Columns whose nodes are all linked are unchanged. TAU <= 0 links every edge present in the
graph. Rationale (Part A): a false homology pair (SPFP) costs FastTree ~5x more than a missed pair (SPFN),
so trading a likely-false merge for an over-split should help trees even if SP error rises.
"""
import argparse
import json
import os

import numpy as np


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def find(par, x):
    root = x
    while par[root] != root:
        root = par[root]
    while par[x] != root:
        par[x], x = root, par[x]
    return root


def split(votedir, tau, out=None, post_field="post_bb"):
    aln = read_fasta(os.path.join(votedir, "out.fasta"))
    E = np.load(os.path.join(votedir, "edges.npz"))
    order = json.load(open(os.path.join(votedir, "subsets.json")))
    # node ids: offset of subset + subset column (MAGUS's order = subsets.json)
    off, o, res_nodes = [], 0, {}
    node_final = []
    for path in order:
        sa = read_fasta(path)
        L = len(next(iter(sa.values())))
        col = np.full(L, -1, dtype=np.int64)
        for t, s in sa.items():
            sc = np.array([i for i, ch in enumerate(s) if ch not in "-."], dtype=np.int64)
            fc = np.array([i for i, ch in enumerate(aln[t]) if ch not in "-."], dtype=np.int64)
            col[sc] = fc
            res_nodes[t] = (o + sc, fc)
        node_final.append(col)
        off.append(o)
        o += L
    node_final = np.concatenate(node_final)
    n = len(node_final)
    a, b, p = E["a"], E["b"], E[post_field]
    same = (node_final[a] == node_final[b]) & (node_final[a] >= 0) & (p > tau)
    par = np.arange(n)
    for x, y in zip(a[same], b[same]):
        rx, ry = find(par, x), find(par, y)
        if rx != ry:
            par[max(rx, ry)] = min(rx, ry)
    comp = np.array([find(par, x) for x in range(n)])
    # new column index: final column major, component (smallest node id) minor
    used = node_final >= 0
    keys = np.unique(np.stack([node_final[used], comp[used]], axis=1), axis=0)
    newcol = {(int(c), int(k)): i for i, (c, k) in enumerate(keys)}
    Lnew = len(keys)
    L = len(next(iter(aln.values())))
    rows = {}
    for t, (nodes, fc) in res_nodes.items():
        row = np.full(Lnew, ord("-"), dtype=np.uint8)
        s = np.frombuffer(aln[t].encode(), dtype=np.uint8)
        idx = np.array([newcol[(int(c), int(comp[x]))] for x, c in zip(nodes, fc)], dtype=np.int64)
        row[idx] = s[fc]
        rows[t] = row.tobytes().decode()
    out = out or os.path.join(votedir, "out.split{}.fasta".format(tau))
    with open(out, "w") as f:
        for t in aln:
            f.write(">{}\n{}\n".format(t, rows[t]))
    return out, L, Lnew


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("votedir")
    ap.add_argument("tau", type=float)
    ap.add_argument("--out")
    ap.add_argument("--post", default="post_bb")
    a = ap.parse_args()
    out, L, Lnew = split(a.votedir, a.tau, a.out, a.post)
    print(json.dumps({"out": out, "cols": L, "cols_split": Lnew}))


if __name__ == "__main__":
    main()
