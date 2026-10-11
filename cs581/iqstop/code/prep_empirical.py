"""Empirical test alignments: random taxon subsamples of curated reference alignments
(Gutell CRW 16S rRNA: 16S.3, 16S.T, 16S.B.ALL; BAliBASE RV100 protein), all-gap columns removed.
Duplicate sequences are left in (IQ-TREE removes them itself)."""
import random, os, sys
D = '/opt/data/Datasets'; out = sys.argv[1]; os.makedirs(out, exist_ok=True)
def read_fa(p):
    seqs, name = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): name = l[1:].split()[0]; seqs[name] = []
        elif name: seqs[name].append(l)
    return {k: ''.join(v).upper() for k, v in seqs.items()}
def write(name, seqs, keep, gapchars='-.?'):
    names = sorted(keep); L = len(seqs[names[0]])
    cols = [j for j in range(L) if any(seqs[n][j] not in gapchars for n in names)]
    with open(f'{out}/{name}.fa', 'w') as f:
        for n in names: f.write(f'>{n}\n' + ''.join(seqs[n][j] for j in cols).replace('.', '-') + '\n')
    print(name, len(names), len(cols))
random.seed(581)
for src in ['16S.3', '16S.T', '16S.B.ALL']:
    s = read_fa(f'{D}/Gutell/{src}/R0/true_align_clean.txt')
    # drop sequences that are mostly gaps (fragments)
    full = [n for n, x in s.items() if sum(c not in '-.?' for c in x) > 1000]
    for k, n in enumerate([100, 150, 200, 300]):
        write(f'emp16S_{src.replace(".","")}_{k}_n{n}', s, random.sample(full, n))
for d in sorted(os.listdir(f'{D}/balibase')):
    s = read_fa(f'{D}/balibase/{d}/model/true.fasta')
    n = min(len(s), random.choice([100, 150, 200]))
    write(f'empAA_{d}_n{n}', s, random.sample(sorted(s), n))
