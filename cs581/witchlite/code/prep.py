#!/usr/bin/env python3
"""Build a backbone/query instance from a reference (true) alignment.

usage: prep.py TRUE_ALN OUTDIR N_BACKBONE N_QUERY READ_LEN SEED

Writes OUTDIR/backbone.fasta (true alignment restricted to the backbone, all-gap columns
removed), backbone.unaln.fasta, queries.fasta (one random window of READ_LEN bp per query
sequence; READ_LEN<=0 = full length), truth.json (for each query residue: its backbone column,
or -1 if its true column holds no backbone residue) and colcount.json (backbone residues per
backbone column). Backbone sequences are drawn from those within 25% of the median length
(as UPP/WITCH do); query sequences from the remaining ones.
"""
import sys, json, random


def read_fasta(path):
    names, seqs, cur = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur:
                    seqs.append(''.join(cur)); cur = []
                names.append(line[1:].split()[0])
            else:
                cur.append(line)
    if cur:
        seqs.append(''.join(cur))
    return names, seqs


def clean(s):
    s = s.upper().replace('U', 'T').replace('.', '-')
    return ''.join(c if c in 'ACGT-' else 'N' for c in s)


def main():
    aln, out, nb, nq, rl, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    import os
    os.makedirs(out, exist_ok=True)
    rng = random.Random(seed)
    names, seqs = read_fasta(aln)
    seqs = [clean(s) for s in seqs]
    lens = [len(s) - s.count('-') for s in seqs]
    med = sorted(lens)[len(lens) // 2]
    ok = [i for i in range(len(names)) if abs(lens[i] - med) <= 0.25 * med]
    rng.shuffle(ok)
    bb = sorted(ok[:nb])
    bbset = set(bb)
    rest = [i for i in range(len(names)) if i not in bbset and lens[i] >= max(rl, 50)]
    rng.shuffle(rest)
    qs = rest[:nq]
    L = len(seqs[0])
    keep = [c for c in range(L) if any(seqs[i][c] != '-' for i in bb)]
    col2bb = {c: j for j, c in enumerate(keep)}
    colcount = [sum(1 for i in bb if seqs[i][c] != '-') for c in keep]
    with open(f'{out}/backbone.fasta', 'w') as f:
        for i in bb:
            f.write(f'>{names[i]}\n' + ''.join(seqs[i][c] for c in keep) + '\n')
    with open(f'{out}/backbone.unaln.fasta', 'w') as f:
        for i in bb:
            f.write(f'>{names[i]}\n' + seqs[i].replace('-', '') + '\n')
    truth = {}
    with open(f'{out}/queries.fasta', 'w') as f:
        for i in qs:
            s = seqs[i]
            cols = [c for c in range(L) if s[c] != '-']
            if rl > 0 and len(cols) > rl:
                st = rng.randrange(0, len(cols) - rl + 1)
                cols = cols[st:st + rl]
            f.write(f'>{names[i]}\n' + ''.join(s[c] for c in cols) + '\n')
            truth[names[i]] = [col2bb.get(c, -1) for c in cols]
    json.dump(truth, open(f'{out}/truth.json', 'w'))
    json.dump(colcount, open(f'{out}/colcount.json', 'w'))
    print(f'backbone {len(bb)} seqs x {len(keep)} cols; {len(qs)} queries; median len {med}')


if __name__ == '__main__':
    main()
