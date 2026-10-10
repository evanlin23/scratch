"""Random k-leaf sub-backbone of a prepared dataset (same queries, alignment columns, model and
true tree), for whole-tree EPA-ng runs that fit in memory.
Usage: python make_subset_dataset.py <datadir> <k> <outdir> [seed=1]"""
import sys, os, random, shutil
import treeswift
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prep_scampp import read_fasta, write_fasta
d, k, o = sys.argv[1], int(sys.argv[2]), sys.argv[3]; seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
os.makedirs(o, exist_ok=True)
t = treeswift.read_tree_newick(f'{d}/rx.raxml.bestTree')
keep = set(random.Random(seed).sample(sorted(v.label for v in t.traverse_leaves()), k))
pt = t.extract_tree_with(keep, suppress_unifurcations=True)
r = pt.root
if len(r.children) == 2:
    a, b = r.children
    m = a if not a.is_leaf() else b
    x = b if m is a else a
    x.edge_length = (x.edge_length or 0) + (m.edge_length or 0)
    r.remove_child(m)
    for c in list(m.children):
        m.remove_child(c); r.add_child(c)
pt.root.edge_length = None; pt.is_rooted = False
pt.write_tree_newick(f'{o}/rx.raxml.bestTree')
s = read_fasta(f'{d}/backbone.fa')
write_fasta(f'{o}/backbone.fa', {n: s[n] for n in sorted(keep)})
for f in ('query_frag.fa', 'query_full.fa', 'queries.txt', 'rx.raxml.bestModel', 'true_tree.txt', 'RAxML_info.REF'):
    if os.path.exists(f'{d}/{f}'):
        shutil.copy(f'{d}/{f}', f'{o}/{f}')
print(o, k)
