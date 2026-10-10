"""Build a batch placement instance from a SCAMPP benchmark dataset (IDB-9257957).

Usage: python prep_scampp.py <dataset_dir> <true_tree> <outdir> [nquery=1000] [seed=1]
dataset_dir holds alignment.fasta, alignment.phylip.raxml.bestTree (backbone topology with
RAxML-ng branch lengths) and alignment.phylip.raxml.bestModel.
Writes (layout of cs581/epang/code/score.py): backbone.fa, query_full.fa, query_frag.fa,
queries.txt, rx.raxml.bestTree (tree pruned to the backbone, unrooted), rx.raxml.bestModel
(site range fixed to the masked length), true_tree.txt.
Queries are removed from the tree all at once (batch placement, as in BSCAMPP), columns with
>95% gaps among backbone sequences are masked, fragments follow prep.py (random start, length
~ N(10% of ungapped length, sd 10 nt)).
"""
import sys, os, random, re
import treeswift
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


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
    dd, true_tree, out = sys.argv[1], os.path.abspath(sys.argv[2]), sys.argv[3]
    nq = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    rng = random.Random(seed)
    os.makedirs(out, exist_ok=True)
    s = read_fasta(f'{dd}/alignment.fasta')
    t = treeswift.read_tree_newick(f'{dd}/alignment.phylip.raxml.bestTree')
    tl = {v.label for v in t.traverse_leaves()}
    names = sorted(set(s) & tl)
    q = set(rng.sample(names, nq))
    bb = [n for n in names if n not in q]
    L = len(s[names[0]])
    keep = [j for j in range(L) if sum(s[n][j] not in '-.?N' for n in bb) >= 0.05 * len(bb)]
    m = {n: ''.join(s[n][j] for j in keep) for n in names}
    write_fasta(f'{out}/backbone.fa', {n: m[n] for n in bb})
    qs = sorted(q)
    write_fasta(f'{out}/query_full.fa', {n: m[n] for n in qs})
    frag = {}
    for n in qs:
        orig = s[n]
        pos = [j for j, c in enumerate(orig) if c != '-']
        ln = min(len(pos), max(20, int(round(rng.gauss(0.1 * len(pos), 10)))))
        st = rng.randint(0, len(pos) - ln)
        lo, hi = pos[st], pos[st + ln - 1]
        fr = ''.join(c if lo <= j <= hi else '-' for j, c in enumerate(orig))
        frag[n] = ''.join(fr[j] for j in keep)
    write_fasta(f'{out}/query_frag.fa', frag)
    open(f'{out}/queries.txt', 'w').write('\n'.join(qs) + '\n')
    # prune the tree to the backbone; keep it unrooted (trifurcating root)
    pt = t.extract_tree_with(set(bb), suppress_unifurcations=True)
    r = pt.root
    if len(r.children) == 2:
        a, b = r.children
        keep_c, merge = (a, b) if not b.is_leaf() else (b, a)
        if merge.is_leaf():
            raise SystemExit('two-leaf tree?')
        # merge child `merge` into root; its edge length moves to keep_c
        keep_c.edge_length = (keep_c.edge_length or 0) + (merge.edge_length or 0)
        r.remove_child(merge)
        for c in list(merge.children):
            merge.remove_child(c); r.add_child(c)
    pt.root.edge_length = None
    pt.write_tree_newick(f'{out}/rx.raxml.bestTree')
    mod = open(f'{dd}/alignment.phylip.raxml.bestModel').read().strip()
    mod = re.sub(r'= *1-\d+', f'= 1-{len(keep)}', mod)
    open(f'{out}/rx.raxml.bestModel', 'w').write(mod + '\n')
    open(f'{out}/true_tree.txt', 'w').write(true_tree + '\n')
    fl = [len(v.replace('-', '')) for v in frag.values()]
    print(f'backbone={len(bb)} queries={nq} sites {L}->{len(keep)} '
          f'pruned_leaves={pt.num_nodes(internal=False)} root_deg={len(pt.root.children)} '
          f'mean_frag_len_after_mask={sum(fl)/len(fl):.1f}')


if __name__ == '__main__':
    main()
