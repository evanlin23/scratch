#!/usr/bin/env python3
"""Second query set for the same held-out genomes: the first L nt after the 341F primer
(CCTACGGGNGGCWGCAG), mimicking the 130-nt single-end reads of the 'indian' dataset.
Usage: loo_queries_v3.py heldout.txt L out.fna"""
import sys, re
held, L, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
iupac = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','M':'[AC]','K':'[GT]','S':'[CG]','W':'[AT]','N':'[ACGT]'}
F = re.compile(''.join(iupac[c] for c in 'CCTACGGGNGGCWGCAG'))
hs = [l.strip() for l in open(held) if l.strip()]
seqs, n = {}, None
for l in open('/opt/mm/root/envs/picrust2/lib/python3.12/site-packages/picrust2/default_files/bacteria/bac_ref/bac_ref.fna'):
    l = l.strip()
    if l.startswith('>'): n = l[1:].split()[0]; seqs[n] = [] if n in hs else None
    elif seqs[n] is not None: seqs[n].append(l)
k = 0
with open(out, 'w') as fo:
    for g in hs:
        u = re.sub(r'[-.]', '', ''.join(seqs[g])).upper()
        m = F.search(u)
        if m and len(u) - m.end() >= L:
            fo.write(f'>{g}\n{u[m.end():m.end() + L]}\n'); k += 1
print(k, 'of', len(hs), 'held-out genomes have 341F', file=sys.stderr)
