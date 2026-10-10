"""Full-data ML baseline for a simulated replicate: IQ-TREE 3 (GTR+G) on the whole alignment.
Usage: python3 sim_fullml.py SIMDIR  -> SIMDIR/iqtree_full.json"""
import json, os, subprocess, sys, time
from phylo import read_tree, fn_fp
d = sys.argv[1]
IQ = os.environ.get("IQTREE", "/opt/mm/root/envs/bio/bin/iqtree3")
if not os.path.exists(f"{d}/iqfull.treefile"):
    t0 = time.time()
    subprocess.run([IQ, "-s", f"{d}/aln.fa", "-m", "GTR+G", "-T", "1", "--prefix", f"{d}/iqfull", "-seed", "1", "-quiet", "-redo"],
                   check=True, capture_output=True)
    sec = time.time() - t0
else:
    sec = None
T = read_tree(f"{d}/true.tre")
r = dict(dir=d, fn_iqtree_full=fn_fp(read_tree(f"{d}/iqfull.treefile"), T)[0], sec=sec)
json.dump(r, open(f"{d}/iqtree_full.json", "w")); print(json.dumps(r), flush=True)
