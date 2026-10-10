"""(cs581/epang/code/prep.py with numpy column masking.) Split an RNASim replicate into backbone + queries, mask gappy sites, make fragments.

Usage: python prep.py <true_align.txt> <outdir> [nquery=1000] [seed=1]
Writes: backbone.fa, query_full.fa, query_frag.fa (all aligned, masked), queries.txt
Fragments follow BSCAMPP: random start, length ~ N(10% of ungapped length, sd 10 nt).
"""
import sys, random, os
import numpy as np

def read_fasta(p):
    seqs, name = {}, None
    for l in open(p):
        l = l.strip()
        if not l:
            continue
        if l.startswith('>'):
            name = l[1:].split()[0]; seqs[name] = []
        else:
            seqs[name].append(l)
    return {k: ''.join(v).upper() for k, v in seqs.items()}

def write_fasta(p, d):
    with open(p, 'w') as f:
        for k, v in d.items():
            f.write(f'>{k}\n{v}\n')

def main():
    aln, out = sys.argv[1], sys.argv[2]
    nq = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    rng = random.Random(seed)
    os.makedirs(out, exist_ok=True)
    s = read_fasta(aln)
    names = sorted(s)
    q = set(rng.sample(names, nq))
    bb = [n for n in names if n not in q]
    L = len(s[names[0]])
    # mask columns with > 95% gaps among backbone sequences (as in BSCAMPP)
    A = np.frombuffer(''.join(s[n] for n in bb).encode(), dtype=np.uint8).reshape(len(bb), L)
    keep = list(np.nonzero((A != ord('-')).sum(0) >= 0.05 * len(bb))[0]); del A
    m = {n: ''.join(s[n][j] for j in keep) for n in names}
    write_fasta(f'{out}/backbone.fa', {n: m[n] for n in bb})
    qs = sorted(q)
    write_fasta(f'{out}/query_full.fa', {n: m[n] for n in qs})
    frag = {}
    for n in qs:
        # fragment on the unmasked sequence, then mask
        orig = s[n]
        pos = [j for j, c in enumerate(orig) if c != '-']
        ln = max(20, int(round(rng.gauss(0.1 * len(pos), 10))))
        st = rng.randint(0, len(pos) - ln)
        lo, hi = pos[st], pos[st + ln - 1]
        fr = ''.join(c if lo <= j <= hi else '-' for j, c in enumerate(orig))
        frag[n] = ''.join(fr[j] for j in keep)
    write_fasta(f'{out}/query_frag.fa', frag)
    open(f'{out}/queries.txt', 'w').write('\n'.join(qs) + '\n')
    fl = [len(v.replace('-', '')) for v in frag.values()]
    print(f'backbone={len(bb)} queries={nq} sites {L}->{len(keep)} mean_frag_len_after_mask={sum(fl)/len(fl):.1f}')

if __name__ == '__main__':
    main()
