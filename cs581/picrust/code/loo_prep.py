#!/usr/bin/env python3
"""Leave-out experiment setup. Pick N reference genomes (seed), cut their 16S V4 amplicon
(515F-Y / 806R primers, primers excluded) from the reference 16S, and write a reduced reference
directory without them (MSA here; tree pruned by loo_prune.R).
Usage: loo_prep.py N seed outdir"""
import sys, re, random, gzip, os, shutil
N, seed, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
R = '/opt/mm/root/envs/picrust2/lib/python3.12/site-packages/picrust2/default_files/bacteria'
iupac = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','M':'[AC]','K':'[GT]','S':'[CG]','W':'[AT]',
         'B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
def rx(p): return ''.join(iupac[c] for c in p)
comp = str.maketrans('ACGTRYMKSWBDHVN', 'TGCAYRKMSWVHDBN')
F = re.compile(rx('GTGYCAGCMGCCGCGGTAA')); Rv = re.compile(rx('GGACTACNVGGGTWTCTAAT'[::-1].translate(comp)))
seqs, order, name = {}, [], None
for l in open(f'{R}/bac_ref/bac_ref.fna'):
    l = l.strip()
    if l.startswith('>'): name = l[1:].split()[0]; seqs[name] = []; order.append(name)
    else: seqs[name].append(l)
seqs = {k: ''.join(v) for k, v in seqs.items()}
amp = {}
for k, s in seqs.items():
    u = re.sub(r'[-.]', '', s).upper()
    f = F.search(u); r = Rv.search(u, f.end() if f else 0) if f else None
    if f and r and 200 <= r.start() - f.end() <= 300:
        amp[k] = u[f.end():r.start()]
print(f'{len(amp)}/{len(seqs)} reference 16S contain both V4 primers', file=sys.stderr)
random.seed(seed)
pick = sorted(random.sample(sorted(amp), N))
os.makedirs(f'{out}/bac_ref', exist_ok=True)
with open(f'{out}/queries.fna', 'w') as fo:
    for k in pick: fo.write(f'>{k}\n{amp[k]}\n')
with open(f'{out}/heldout.txt', 'w') as fo: fo.write('\n'.join(pick) + '\n')
ps = set(pick)
with open(f'{out}/bac_ref/bac_ref.fna', 'w') as fo:
    for k in order:
        if k not in ps: fo.write(f'>{k}\n{seqs[k]}\n')
for ext in ('hmm', 'model', 'raxml_info'):
    shutil.copy(f'{R}/bac_ref/bac_ref.{ext}', f'{out}/bac_ref/bac_ref.{ext}')
print('amplicon length mean', sum(len(amp[k]) for k in pick) / N, file=sys.stderr)
