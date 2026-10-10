#!/usr/bin/env python3
"""Query-to-backbone homology scoring (exact SPFN/SPFP restricted to query-backbone pairs).

Each method is reduced to a map: query residue -> backbone column (or -1 = insertion / unaligned).
Since the backbone alignment is fixed and equals the true alignment restricted to the backbone,
a query residue at true backbone column t has colcount[t] true homologies with backbone residues;
it shares all of them if placed in column t and none otherwise; it predicts colcount[e] pairs
when placed in column e.

usage as a library: load_instance(dir), map_from_alignment(path, backbone_names), metrics(...)
"""
import json, sys
import numpy as np


def read_fasta(path):
    d, name, cur = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    d[name] = ''.join(cur)
                name, cur = line[1:].split()[0], []
            else:
                cur.append(line.strip())
    if name is not None:
        d[name] = ''.join(cur)
    return d


def load_instance(inst):
    truth = json.load(open(f'{inst}/truth.json'))
    colcount = np.array(json.load(open(f'{inst}/colcount.json')))
    bb = read_fasta(f'{inst}/backbone.fasta')
    return truth, colcount, bb


def map_from_alignment(path, bbnames, queries):
    """Map query residues of an extended alignment (backbone + queries) to backbone columns."""
    aln = read_fasta(path)
    bbs = [aln[n] for n in bbnames]
    L = len(bbs[0])
    arr = np.frombuffer(''.join(bbs).encode(), dtype=np.uint8).reshape(len(bbs), L)
    has = ((arr != ord('-')) & (arr != ord('.'))).any(axis=0)
    colmap = np.full(L, -1)
    colmap[has] = np.arange(has.sum())
    out = {}
    for q in queries:
        s = aln.get(q)
        if s is None:
            out[q] = None
            continue
        res = []
        for c, ch in enumerate(s):
            if ch not in '-.':
                res.append(int(colmap[c]) if ch.isupper() else -1)
        out[q] = res
    return out, int(has.sum())


def metrics(truth, colcount, est):
    """Returns pooled (spfn, spfp) and per-query arrays (spfn_q, spfp_q)."""
    names = list(truth.keys())
    S = T = E = 0
    fn, fp = [], []
    for q in names:
        t = np.array(truth[q])
        e = est.get(q)
        if e is None or len(e) != len(t):
            e = np.full(len(t), -1) if e is None else None
            if e is None:
                raise ValueError(f'length mismatch for {q}')
        e = np.asarray(e)
        tc = np.where(t >= 0, colcount[np.maximum(t, 0)], 0)
        ec = np.where(e >= 0, colcount[np.maximum(e, 0)], 0)
        sh = np.where((t == e) & (t >= 0), tc, 0)
        s, tt, ee = sh.sum(), tc.sum(), ec.sum()
        S += s; T += tt; E += ee
        fn.append(1 - s / tt if tt else 0.0)
        fp.append(1 - s / ee if ee else 0.0)
    return 1 - S / T, (1 - S / E if E else 0.0), np.array(fn), np.array(fp)


if __name__ == '__main__':
    inst, aln = sys.argv[1], sys.argv[2]
    truth, cc, bb = load_instance(inst)
    est, nb = map_from_alignment(aln, list(bb.keys()), truth.keys())
    assert nb == len(cc), (nb, len(cc))
    fn, fp, _, _ = metrics(truth, cc, est)
    print(f'SPFN {fn:.4f} SPFP {fp:.4f}')
