"""GTR+G parameters for a large backbone from RAxML-NG --evaluate on a random k-leaf pruned
subtree; writes rx.raxml.bestModel (site range = full alignment) and rx.raxml.bestTree (=
the FastTree tree, unrooted). Usage: python model_from_subtree.py <dir> [k=5000] [threads=1]"""
import sys, os, random, subprocess, re
import treeswift
d = sys.argv[1]; k = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
T = sys.argv[3] if len(sys.argv) > 3 else '1'
P = '/opt/mm/root/envs/place/bin'
t = treeswift.read_tree_newick(f'{d}/fasttree.tre')
t.root.edge_length = None
for v in t.traverse_internal():
    v.label = None  # drop FastTree support values
t.is_rooted = False
t.write_tree_newick(f'{d}/rx.raxml.bestTree')
names = [v.label for v in t.traverse_leaves()]
sub = set(random.Random(3).sample(names, k))
st = t.extract_tree_with(sub, suppress_unifurcations=True)
r = st.root
if len(r.children) == 2:  # unroot
    a, b = r.children
    m = a if not a.is_leaf() else b
    o = b if m is a else a
    o.edge_length = (o.edge_length or 0) + (m.edge_length or 0)
    r.remove_child(m)
    for c in list(m.children):
        m.remove_child(c); r.add_child(c)
st.root.edge_length = None
st.is_rooted = False; t.is_rooted = False
st.write_tree_newick(f'{d}/msub.tre')
L = 0
with open(f'{d}/backbone.fa') as f, open(f'{d}/msub.fa', 'w') as o:
    keep = False
    for l in f:
        if l.startswith('>'):
            keep = l[1:].strip() in sub
        elif not L:
            L = len(l.strip())
        if keep:
            o.write(l)
subprocess.run([f'{P}/raxml-ng', '--evaluate', '--msa', f'{d}/msub.fa', '--tree', f'{d}/msub.tre',
                '--model', 'GTR+G', '--prefix', f'{d}/msub', '--threads', T, '--redo'], check=True)
m = open(f'{d}/msub.raxml.bestModel').read().strip()
open(f'{d}/rx.raxml.bestModel', 'w').write(re.sub(r'= *1-\d+', f'= 1-{L}', m) + '\n')
print(open(f'{d}/rx.raxml.bestModel').read())
