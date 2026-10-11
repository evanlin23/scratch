#!/usr/bin/env python3
"""Placement delta error for query alignments (EPA-ng on the fixed FastTree backbone tree).

For each method's query->backbone-column map, build the query MSA, place with EPA-ng (best LWR
edge), and compute per-query delta error  = FN(T+q, T*|B+q) - FN(T, T*|B)  (missing true
bipartitions, counts), T* = true tree. Methods are given as name=path pairs; path is either a
WITCH output alignment or a JSON map (BLAST).
usage: placement.py INST TRUE_TREE OUTDIR name=path [name=path ...]
"""
import sys, os, re, json, subprocess
import numpy as np
import treeswift
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score import load_instance, map_from_alignment

EPA = '/opt/mm/root/envs/epa/bin/epa-ng'
inst, ttree, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
truth, cc, bb = load_instance(inst)
bbn = list(bb); L = len(cc)
qn = list(truth)
idx = {n: i for i, n in enumerate(bbn + qn)}
BM = sum(1 << idx[n] for n in bbn)


def splits(tree, mask):
    """nontrivial bipartitions (canonical bitmasks within `mask`) of a treeswift tree"""
    out = set()
    nl = bin(mask).count('1')
    for nd in tree.traverse_postorder():
        if nd.is_leaf():
            nd.m = (1 << idx[nd.label]) if nd.label in idx else 0
        else:
            nd.m = 0
            for c in nd.children:
                nd.m |= c.m
        s = nd.m & mask
        k = bin(s).count('1')
        if 2 <= k <= nl - 2:
            if s & 1:
                s = mask & ~s
            out.add(s)
    return out


# true tree restricted to backbone + all queries
T = treeswift.read_tree_newick(ttree)
keep = set(idx)
T = T.extract_tree_with(keep)
true_all = list(splits(T, (1 << len(idx)) - 1))

def true_splits(mask):
    nl = bin(mask).count('1'); out = set()
    for s in true_all:
        s = s & mask
        k = bin(s).count('1')
        if 2 <= k <= nl - 2:
            if s & 1:
                s = mask & ~s
            out.add(s)
    return out

bbtree = treeswift.read_tree_newick(f'{inst}/backbone.tre')
base_fn = len(true_splits(BM) - splits(bbtree, BM))
print('backbone tree FN (count):', base_fn, flush=True)

res = {}
for arg in sys.argv[4:]:
    name, path = arg.split('=', 1)
    od = f'{out}/{name}'
    os.makedirs(od, exist_ok=True)
    est = json.load(open(path)) if path.endswith('.json') else map_from_alignment(path, bbn, truth)[0]
    seqs = {}
    cur = None
    for line in open(f'{inst}/queries.fasta'):
        line = line.strip()
        if line.startswith('>'):
            cur = line[1:].split()[0]; seqs[cur] = ''
        else:
            seqs[cur] += line
    with open(f'{od}/q.fa', 'w') as f:
        for q in qn:
            row = ['-'] * L
            for r, c in enumerate(est.get(q) or []):
                if c >= 0:
                    row[c] = seqs[q][r]
            if any(ch != '-' for ch in row):  # EPA-ng aborts on all-gap queries: those count as unplaced
                f.write(f'>{q}\n' + ''.join(row) + '\n')
    if not os.path.exists(f'{od}/epa_result.jplace'):
        subprocess.run([EPA, '--ref-msa', f'{inst}/backbone.fasta', '--tree', f'{inst}/backbone.tre',
                        '--query', f'{od}/q.fa', '--model', 'GTR+G', '-T', '4', '-w', od, '--redo'],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    jp = json.load(open(f'{od}/epa_result.jplace'))
    fields = jp['fields']
    tstr = re.sub(r'([^(),:;]*):([0-9.eE+-]+)\{(\d+)\}', lambda m: f'{m.group(1)}@{m.group(3)}:{m.group(2)}', jp['tree'])
    deltas = []
    best = {}
    for p in jp['placements']:
        names = [n[0] if isinstance(n, list) else n for n in p.get('n', p.get('nm', []))]
        row = max(p['p'], key=lambda r: r[fields.index('like_weight_ratio')])
        for nm in names:
            best[nm] = row[fields.index('edge_num')]
    for q in qn:
        if q not in best:
            deltas.append(np.nan); continue
        t = treeswift.read_tree_newick(tstr)
        tgt = None
        for nd in t.traverse_preorder():
            lab = nd.label or ''
            if '@' in lab:
                l, e = lab.rsplit('@', 1)
                nd.label = l if l else None
                if int(e) == best[q]:
                    tgt = nd
        # insert q on the edge above tgt
        par = tgt.parent
        par.remove_child(tgt)
        mid = treeswift.Node(); par.add_child(mid); mid.add_child(tgt)
        mid.add_child(treeswift.Node(label=q))
        M = BM | (1 << idx[q])
        fn = len(true_splits(M) - splits(t, M))
        deltas.append(fn - base_fn)
    d = np.array(deltas, dtype=float)
    res[name] = dict(mean_delta=float(np.nanmean(d)), unplaced=int(np.isnan(d).sum()), delta_q=d.tolist())
    print(f'{name:20s} mean delta error {np.nanmean(d):.3f} (unplaced {int(np.isnan(d).sum())})', flush=True)
    json.dump(res, open(f'{out}/placement.json', 'w'))
