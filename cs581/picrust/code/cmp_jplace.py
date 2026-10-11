#!/usr/bin/env python3
"""Compare two EPA-ng jplace files (same queries, same reference tree).
Per query: best edge (max LWR), LWR, edge distance between best edges (number of edges on the
tree path), patristic distance between placement points, # placements kept.
Usage: cmp_jplace.py A.jplace B.jplace out.tsv"""
import json, re, sys, collections

def parse_tree(s):
    # returns parent map over edge numbers and edge lengths; nodes are identified by edge_num of
    # the edge above them; root gets -1
    parent, blen, stack, cur = {}, {}, [], []
    tok = re.finditer(r'\(|\)|,|;|[^(),;]+', s)
    children = collections.defaultdict(list)
    pending = []  # stack of child lists
    last = None
    for m in tok:
        t = m.group()
        if t == '(':
            pending.append([])
        elif t == ',':
            pending[-1].append(last); last = None
        elif t == ')':
            pending[-1].append(last); kids = pending.pop(); last = ('inner', kids)
        elif t == ';':
            break
        else:
            # label[:len]{edge}
            mm = re.match(r'([^:{]*)(?::([^{]*))?(?:\{(\d+)\})?', t.strip())
            lab, l, e = mm.groups()
            if last is None or last[0] != 'inner':
                last = ('leaf', [])
            node = {'kids': last[1], 'label': lab, 'len': float(l) if l else 0.0,
                    'edge': int(e) if e is not None else -1}
            last = ('node', node)
            # fix kids for leaf/inner
    root = last[1] if last[0] == 'node' else {'kids': last[1], 'label': '', 'len': 0, 'edge': -1}
    # build parent map by edge num
    par, ln, depth, name = {}, {}, {}, {}
    def walk(n, p, d):
        e = n['edge']; par[e] = p; ln[e] = n['len']; depth[e] = d; name[e] = n['label']
        for k in n['kids']:
            walk(k[1], e, d + 1)
    sys.setrecursionlimit(1000000)
    walk(root, None, 0)
    return par, ln, depth, name

def path(par, depth, ln, a, b, pa, pb):
    """edges on the path between edge a and edge b (count) and patristic dist between points at
    distal length pa (from the lower node of a, toward its parent) and pb."""
    if a == b:
        return 0, abs(pa - pb)
    # distance from point on edge a up to the top node of a: ln[a]-pa ; from point to lower node: pa
    x, y, n = a, b, 0
    dx, dy = ln[a] - pa, ln[b] - pb  # distances from points to the upper ends
    ex, ey = 0, 0
    while depth[x] > depth[y]:
        x = par[x]; ex += 1
        if x != a and par[x] is not None and x != -1: pass
    while depth[y] > depth[x]:
        y = par[y]; ey += 1
    while x != y:
        x = par[x]; y = par[y]; ex += 1; ey += 1
    # x==y is lowest common ancestor node (identified by its edge). Count edges strictly between.
    if x == b:  # a is below b
        return ex, None
    if x == a:
        return ey, None
    return ex + ey - 1, None

def best(d):
    out = {}
    f = d['fields']; ie, il, iw = f.index('edge_num'), f.index('likelihood'), f.index('like_weight_ratio')
    idl, ipl = f.index('distal_length'), f.index('pendant_length')
    for p in d['placements']:
        b = max(p['p'], key=lambda r: r[iw])
        names = p.get('n') or [x[0] for x in p['nm']]
        for n in names:
            out[n] = dict(edge=b[ie], logl=b[il], lwr=b[iw], distal=b[idl], pendant=b[ipl], npl=len(p['p']))
    return out

A, B = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
par, ln, depth, _ = parse_tree(A['tree'])
a, b = best(A), best(B)
with open(sys.argv[3], 'w') as fo:
    fo.write('query\tedgeA\tedgeB\tedge_dist\tlwrA\tlwrB\tloglA\tloglB\tpendA\tpendB\tnplA\tnplB\n')
    for q in sorted(a):
        x, y = a[q], b[q]
        ed, _ = path(par, depth, ln, x['edge'], y['edge'], x['distal'], y['distal'])
        fo.write(f"{q}\t{x['edge']}\t{y['edge']}\t{ed}\t{x['lwr']:.4f}\t{y['lwr']:.4f}\t{x['logl']:.3f}\t{y['logl']:.3f}\t{x['pendant']:.5f}\t{y['pendant']:.5f}\t{x['npl']}\t{y['npl']}\n")
