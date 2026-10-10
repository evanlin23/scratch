"""Offline replay of stopping rules on IQ-TREE 3 per-iteration traces.
For each dataset: parse iqdef1.itertrace, collect every new-best tree + comparators, evaluate all of them
with one IQ-TREE -z run (same model re-estimated on tree 1, branch lengths re-optimized per tree, site
log-likelihoods + AU test), then replay rules and write results/<ds>.json.
Usage: python replay.py <dataset> [<dataset> ...]"""
import os, sys, json, math, subprocess, re
import numpy as np
from scipy.stats import norm
IQ = '/opt/work/bin/iqtree3'; W = '/opt/work/runs'
OUT = '/home/user/scratch/cs581/iqstop/results/per_dataset'

def aln_orig(ds):
    return f'/opt/work/sim/{ds}.phy' if ds.startswith('rg') else f'/opt/work/emp/{ds}.fa'
def aln_of(ds):
    u = f'{W}/{ds}/iqdef1.uniqueseq.phy'
    return u if os.path.exists(u) else aln_orig(ds)

def taxa_of(path):
    if path.endswith('.phy'):
        return [l.split()[0] for l in open(path).read().strip().split('\n')[1:] if l.strip()]
    return [l[1:].split()[0] for l in open(path) if l.startswith('>')]

def relabel(tr, names):
    return re.sub(r'([(,])(\d+):', lambda m: m.group(1) + names[int(m.group(2))] + ':', tr)

def prune(tr, keep):
    import dendropy
    t = dendropy.Tree.get(data=tr, schema='newick', preserve_underscores=True)
    drop = [n.taxon.label for n in t.leaf_node_iter() if n.taxon.label not in keep]
    if drop: t.prune_taxa_with_labels(drop)
    return t.as_string(schema='newick', suppress_rooting=True).strip()

def model_of(ds):
    return 'LG+G4' if ds.startswith('empAA') else 'GTR+F+I+G4'

def read_trace(p):
    rows = []
    for l in open(p):
        it, t, sc, best, newtopo, tr = l.rstrip('\n').split('\t')
        # keep the tree only for new-best *topologies* (IQ-TREE resets its -nstop counter only on these)
        rows.append((int(it), float(t), float(sc), float(best), tr if newtopo == '1' else None))
    return rows

def wall(ds, tag):
    p = f'{W}/{ds}/{tag}.done'
    if not os.path.exists(p): return None
    f = open(p).read().split(); return float(f[3]) if f[2] == '0' else None

def evaluate(ds, trees):
    d = f'{W}/{ds}/eval'; os.makedirs(d, exist_ok=True)
    content = ''.join(t.strip() + '\n' for t in trees)
    same = os.path.exists(f'{d}/all.nwk') and open(f'{d}/all.nwk').read() == content
    open(f'{d}/all.nwk', 'w').write(content)
    if not (same and os.path.exists(f'{d}/ev.sitelh')):
        subprocess.run([IQ, '-s', aln_of(ds), '-m', model_of(ds), '-z', f'{d}/all.nwk', '-n', '0', '-wsl',
                        '-zb', '1000', '-au', '-T', '1', '-seed', '1', '-pre', f'{d}/ev', '-redo', '--quiet'], check=True)
    L = []
    for l in open(f'{d}/ev.sitelh').read().split('\n')[1:]:
        if l.strip(): L.append([float(x) for x in l.split()[1:]])
    L = np.array(L)
    # AU p-values from .iqtree table
    au = {}
    txt = open(f'{d}/ev.iqtree').read()
    sec = txt[txt.find('USER TREES'):]
    for l in sec.split('\n'):
        m = re.match(r'\s*(\d+)\s+(-?[\d.]+)\s+([\d.]+)\s+(.*)', l)
        if m:
            toks = m.group(4).replace('+', ' ').replace('-', ' ').split()
            au[int(m.group(1)) - 1] = float(toks[-1]) if toks else None
    return L, au

def kh_better(sl_new, sl_ref, alpha, eps=0.0):
    """one-sided fast KH (normal approx, Goldman et al. 2000): is new significantly better than ref?"""
    d = sl_new - sl_ref; n = len(d); delta = d.sum(); sd = math.sqrt(n * d.var(ddof=1)) if n > 1 else 0
    return delta > max(norm.ppf(1 - alpha) * sd, eps), delta, sd

def replay(rows, idx_of_it, L, rule):
    """returns stop iteration. rows: trace; idx_of_it: iteration -> index (row in L) of best tree at that iteration."""
    kind = rule['kind']; last_it = rows[-1][0]
    its = [r[0] for r in rows]
    impr = [r[0] for r in rows if r[4] is not None]
    if kind == 'default':
        return last_it
    if kind == 'nstop':
        N = rule['N']; last = 0
        for it in its:
            if it in impr: last = it
            if it > last + N: return it
        return last_it
    if kind == 'khpatience':   # count only KH-significant improvements as "success"
        N, a = rule['N'], rule['alpha']; last = 0; ref = None
        for it in its:
            if it in impr:
                k = idx_of_it[it]
                if ref is None or kh_better(L[:, k], L[:, ref], a, rule.get('eps', 0))[0]:
                    last, ref = it, k
            if it > last + N: return it
        return last_it
    if kind == 'khwindow':     # every k iterations after warm-up, test best now vs best k its ago; stop if not significant
        k, a, warm = rule['k'], rule['alpha'], rule.get('warm', 20)
        prev = idx_of_it[warm] if warm in idx_of_it else None
        for it in its:
            if it > warm and (it - warm) % k == 0:
                cur = idx_of_it[it]
                if rule.get('mult'):
                    m = sum(1 for j in impr if it - k < j <= it)
                    aa = a / max(m, 1)
                else: aa = a
                if cur == prev or not kh_better(L[:, cur], L[:, prev], aa, rule.get('eps', 0))[0]:
                    return it
                prev = cur
        return last_it
    raise ValueError(kind)

RULES = {
    'default(nstop100)': dict(kind='default'),
    'nstop50': dict(kind='nstop', N=50), 'nstop20': dict(kind='nstop', N=20), 'nstop10': dict(kind='nstop', N=10),
    'KHpat100': dict(kind='khpatience', N=100, alpha=0.05),
    'KHpat50': dict(kind='khpatience', N=50, alpha=0.05),
    'KHpat20': dict(kind='khpatience', N=20, alpha=0.05),
    'KHwin10': dict(kind='khwindow', k=10, alpha=0.05),
    'KHwin20': dict(kind='khwindow', k=20, alpha=0.05),
    'KHwin10mult': dict(kind='khwindow', k=10, alpha=0.05, mult=True),
}

def main(ds):
    names = taxa_of(aln_of(ds)); keep = set(names)
    rows = [(a, b, c, d, relabel(e, names) if e else None) for a, b, c, d, e in read_trace(f'{W}/{ds}/iqdef1.itertrace')]
    T = wall(ds, 'iqdef1'); post = T - rows[-1][1]
    trees, tid = [], {}
    def add(t):
        if t not in tid: tid[t] = len(trees); trees.append(t)
        return tid[t]
    idx_of_it, cur = {}, None
    for it, t, sc, best, tr in rows:
        if tr is not None: cur = add(tr)
        idx_of_it[it] = cur
    comps = {}
    for tag, fn in [('iqdef1_final', 'iqdef1.treefile'), ('iqdef2', 'iqdef2.treefile'), ('iqfast', 'iqfast.treefile'),
                    ('rxfast', 'rxfast.raxml.bestTree')]:
        p = f'{W}/{ds}/{fn}'
        if os.path.exists(p): comps[tag] = add(prune(open(p).read().strip(), keep))
    L, au = evaluate(ds, trees)
    lnl = L.sum(0)
    res = dict(ds=ds, ntaxa=None, nsites=L.shape[0], wall_default=T, post=post, last_it=rows[-1][0],
               n_impr=sum(1 for r in rows if r[4] is not None), rules={}, comps={}, trees=trees)
    for name, rule in RULES.items():
        s = replay(rows, idx_of_it, L, rule)
        t_s = [r[1] for r in rows if r[0] == s][0]
        k = idx_of_it[s]
        res['rules'][name] = dict(stop_it=s, time=t_s + post, tree=k, lnl=float(lnl[k]), au=au.get(k))
    for tag, k in comps.items():
        res['comps'][tag] = dict(tree=k, lnl=float(lnl[k]), au=au.get(k), time=wall(ds, tag.replace('_final', '')))
    os.makedirs(OUT, exist_ok=True)
    json.dump(res, open(f'{OUT}/{ds}.json', 'w'))
    print(ds, 'done', {n: (r['stop_it'], round(r['time'], 1), round(r['lnl'] - res['rules']['default(nstop100)']['lnl'], 2)) for n, r in res['rules'].items()})

if __name__ == '__main__':
    for ds in sys.argv[1:]:
        try: main(ds)
        except Exception as e: print(ds, 'FAILED', repr(e))
