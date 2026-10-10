"""Simulate DNA alignments with AliSim on RAxML Grove empirical trees + their fitted GTR+(I)+G params.
Usage: python simulate.py N outdir   (templates stratified by #taxa from /tmp/rg_cand.txt)"""
import os, re, sys, random, subprocess
RG = '/tmp/rg/trees'; IQ = '/opt/mm/root/envs/bio/bin/iqtree3'
N, out = int(sys.argv[1]), sys.argv[2]
cand = [l.split() for l in open('/tmp/rg_cand.txt')]
cand = [(d, int(t), int(s)) for d, t, s in cand]
random.seed(581)
bins = [(50, 100), (100, 200), (200, 500)]
picked = []
for lo, hi in bins:
    c = [x for x in cand if lo <= x[1] < hi]
    picked += random.sample(c, N // len(bins))
os.makedirs(out, exist_ok=True)
for d, nt, ns in picked:
    s = open(f'{RG}/{d}/log_0.txt').read()
    rates = [float(x) for x in re.findall(r'rate [ACGT] <-> [ACGT]: ([\d.]+)', s)][:6]
    freqs = [float(x) for x in re.findall(r'freq pi\([ACGT]\): ([\d.]+)', s)][:4]
    alpha = float(re.search(r'alpha: ([\d.]+)', s).group(1))
    inv = re.search(r'invar: ([\d.]+)', s)
    m = 'GTR{%s}+F{%s}' % (','.join(map(str, rates[:5])), ','.join(map(str, freqs)))
    if inv: m += '+I{%s}' % inv.group(1)
    m += '+G4{%s}' % alpha
    name = f'rg{d}_n{nt}'
    if os.path.exists(f'{out}/{name}.phy'): continue
    # strip zero-length polytomy issues: AliSim accepts the RAxML tree as is
    subprocess.run([IQ, '--alisim', f'{out}/{name}', '-t', f'{RG}/{d}/tree_best.newick', '-m', m,
                    '--length', str(ns), '-seed', '1', '-af', 'phy', '-redo', '--quiet'], check=True)
    open(f'{out}/{name}.model', 'w').write(m + '\n')
    os.system(f'cp {RG}/{d}/tree_best.newick {out}/{name}.true.nwk')
    print(name, nt, ns, m)
