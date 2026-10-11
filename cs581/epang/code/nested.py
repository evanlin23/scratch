"""Nested-subtree experiment: place the same queries into subtrees of growing size.

Step 'prep': choose centers, nested subtrees (BSCAMPP's own edge-length expansion around
the center leaf), and the queries for each center (queries whose closest backbone leaf by
Hamming distance on the true alignment lies among the INNER leaves nearest the center).
Writes <d>/nested/c<i>/k<k>/{tree.nwk,ref.fa,full.fa,frag.fa}.

Usage: python nested.py prep <datadir> [ncenters=20] [inner=200] [sizes=500,1000,2000,3000,5000]
"""
import sys, os, random, json
import numpy as np
import treeswift
from bscampp import utils as bu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prep import read_fasta, write_fasta


def closest_leaf(ref_names, R, qseq):
    """closest reference leaf by Hamming distance over the query's non-gap sites"""
    q = np.frombuffer(qseq.encode(), dtype=np.uint8)
    mask = q != ord('-')
    d = (R[:, mask] != q[mask]).sum(1)
    return ref_names[int(d.argmin())]


def prep(d, ncent=20, inner=200, sizes=(500, 1000, 2000, 3000, 5000)):
    rng = random.Random(7)
    ref = read_fasta(f'{d}/backbone.fa')
    full = read_fasta(f'{d}/query_full.fa')
    frag = read_fasta(f'{d}/query_frag.fa')
    names = sorted(ref)
    R = np.array([np.frombuffer(ref[n].encode(), dtype=np.uint8) for n in names])
    tree = treeswift.read_tree_newick(f'{d}/rx.raxml.bestTree')
    leaf = tree.label_to_node(selection='leaves')
    # seed leaf of each query (from the fragment, as a real pipeline would see it)
    seed = {q: closest_leaf(names, R, frag[q]) for q in sorted(frag)}
    centers = rng.sample(names, ncent)
    used = set()
    meta = {'centers': {}, 'sizes': list(sizes)}
    for i, c in enumerate(centers):
        inner_set = set(bu.subtree_nodes_with_edge_length(tree, leaf[c], inner))
        qs = [q for q in sorted(seed) if seed[q] in inner_set and q not in used]
        used.update(qs)
        if not qs:
            continue
        meta['centers'][f'c{i}'] = {'center': c, 'queries': qs}
        for k in sizes:
            labels = bu.subtree_nodes_with_edge_length(tree, leaf[c], k)
            sub = tree.extract_tree_with(set(labels))
            sub.resolve_polytomies(); sub.suppress_unifurcations()
            od = f'{d}/nested/c{i}/k{k}'
            os.makedirs(od, exist_ok=True)
            sub.write_tree_newick(f'{od}/tree.nwk', hide_rooted_prefix=True)
            write_fasta(f'{od}/ref.fa', {n: ref[n] for n in labels})
            write_fasta(f'{od}/full.fa', {q: full[q] for q in qs})
            write_fasta(f'{od}/frag.fa', {q: frag[q] for q in qs})
    allq = sorted(used)
    meta['queries'] = allq
    os.makedirs(f'{d}/nested/all', exist_ok=True)
    write_fasta(f'{d}/nested/all/full.fa', {q: full[q] for q in allq})
    write_fasta(f'{d}/nested/all/frag.fa', {q: frag[q] for q in allq})
    json.dump(meta, open(f'{d}/nested/meta.json', 'w'), indent=1)
    json.dump(seed, open(f'{d}/seed.json', 'w'))
    print(f'{len(meta["centers"])} centers, {len(allq)} queries')


if __name__ == '__main__':
    if sys.argv[1] == 'prep':
        a = sys.argv[2:]
        prep(a[0], *(int(x) for x in a[1:3]),
             **({'sizes': tuple(int(x) for x in a[3].split(','))} if len(a) > 3 else {}))
