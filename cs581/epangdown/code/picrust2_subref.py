"""PICRUSt2 bacterial reference (26,868 tips; EPA-ng needs >8 GB, est. ~40 GB, so it does not
fit on this 15 GB machine) reduced to a random k-tip subset: writes <out>/sub_ref.{tre,fna,model,hmm}.
Usage: python picrust2_subref.py <k> <outdir> [seed=1]"""
import sys, os, random, shutil
import treeswift
R = '/opt/mm/root/envs/picrust2/lib/python3.11/site-packages/picrust2/default_files/bacteria/bac_ref'
k, o = int(sys.argv[1]), sys.argv[2]; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
os.makedirs(o, exist_ok=True)
t = treeswift.read_tree_newick(f'{R}/bac_ref.tre')
names = sorted(v.label for v in t.traverse_leaves())
keep = set(random.Random(seed).sample(names, k))
st = t.extract_tree_with(keep, suppress_unifurcations=True)
for v in st.traverse_internal():
    v.label = None
r = st.root
if len(r.children) == 2:
    a, b = r.children
    m = a if not a.is_leaf() else b
    x = b if m is a else a
    x.edge_length = (x.edge_length or 0) + (m.edge_length or 0)
    r.remove_child(m)
    for c in list(m.children):
        m.remove_child(c); r.add_child(c)
st.root.edge_length = None; st.is_rooted = False
st.write_tree_newick(f'{o}/sub_ref.tre')
with open(f'{R}/bac_ref.fna') as f, open(f'{o}/sub_ref.fna', 'w') as g:
    w = False
    for l in f:
        if l.startswith('>'):
            w = l[1:].strip() in keep
        if w:
            g.write(l)
shutil.copy(f'{R}/bac_ref.model', f'{o}/sub_ref.model')
# HMM rebuilt from the subset (hmmalign --mapali checks the MSA checksum)
import subprocess
subprocess.run(["/opt/mm/root/envs/picrust2/bin/hmmbuild", "--dna", f"{o}/sub_ref.hmm", f"{o}/sub_ref.fna"], check=True, stdout=subprocess.DEVNULL)
# default bacterial trait tables restricted to the kept tips
import gzip
for tname in ('16S', 'ko', 'ec'):
    with gzip.open(f'{R}/../{tname}.txt.gz', 'rt') as f, gzip.open(f'{o}/{tname}.txt.gz', 'wt') as g:
        for i, l in enumerate(f):
            if i == 0 or l.split('\t', 1)[0] in keep:
                g.write(l)
print(k, 'tips')
