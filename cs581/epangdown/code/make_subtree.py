"""Make a k-leaf subtree instance from a prepared dataset dir: the k leaves that come first in a
preorder traversal of the backbone tree (a union of whole clades), pruned and unrooted.
Usage: python make_subtree.py <datadir> <k> <outdir> <nq>   (writes tree.nwk, ref.fa, q_frag.fa, q_full.fa)"""
import sys, os, random
import treeswift
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prep_scampp import read_fasta, write_fasta
d, k, o, nq = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
os.makedirs(o, exist_ok=True)
t = treeswift.read_tree_newick(f'{d}/rx.raxml.bestTree')
leaves = [v.label for v in t.traverse_preorder() if v.is_leaf()][:k]
pt = t.extract_tree_with(set(leaves), suppress_unifurcations=True)
r = pt.root
if len(r.children) == 2:
    a, b = r.children
    m = a if not a.is_leaf() else b
    other = b if m is a else a
    other.edge_length = (other.edge_length or 0) + (m.edge_length or 0)
    r.remove_child(m)
    for c in list(m.children):
        m.remove_child(c); r.add_child(c)
pt.root.edge_length = None; pt.is_rooted = False
pt.write_tree_newick(f'{o}/tree.nwk')
s = read_fasta(f'{d}/backbone.fa')
write_fasta(f'{o}/ref.fa', {n: s[n] for n in leaves})
rng = random.Random(7)
for qt in ('frag', 'full'):
    qq = read_fasta(f'{d}/query_{qt}.fa')
    keys = sorted(qq)[:nq]
    write_fasta(f'{o}/q_{qt}.fa', {n: qq[n] for n in keys})
print(o, len(leaves), 'leaves')
