"""TreeMerge (Molloy & Warnow 2019; original github.com/ekmolloy/treemerge with PAUP*) on a
simulated replicate, using topological (node-count) distances on the guide tree, as in
Park et al. 2021. Usage: python3 sim_treemerge.py SIMDIR STMODE -> SIMDIR/sub_STMODE/treemerge.json
Needs: Python 2.7 env with dendropy 4.3 + networkx (TM_PY), PAUP* (PAUP), TreeMerge (TM_DIR)."""
import glob, json, os, shutil, subprocess, sys, time
from phylo import read_tree, fn_fp
TM_PY = os.environ.get("TM_PY", "/opt/mm/root/envs/tm27/bin/python")
TM_DIR = os.environ.get("TM_DIR", "/opt/tools/treemerge/python")
PAUP = os.environ.get("PAUP", "/opt/tools/paup/paup")
d, st = sys.argv[1], sys.argv[2]
sd = f"{d}/sub_{st}"
G = read_tree(f"{d}/guide.tre")
names = sorted(G.label.values())
leafid = {G.label[v]: v for v in G.label}
# node (edge-count) distances on the guide tree
D = {}
for a in names:
    dist = {leafid[a]: 0}
    st_ = [leafid[a]]
    while st_:
        v = st_.pop()
        for u in G.adj[v]:
            if u not in dist:
                dist[u] = dist[v] + 1
                st_.append(u)
    D[a] = [dist[leafid[b]] for b in names]
mat = f"{sd}/guide_node.mat"
with open(mat, "w") as f:
    f.write(f"{len(names)}\n")
    for a in names:
        f.write(a + " " + " ".join(map(str, D[a])) + "\n")
open(mat + "_taxlist", "w").write("\n".join(names) + "\n")
work = f"{sd}/tm_work"
shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
out = f"{sd}/treemerge.tre"
if os.path.exists(out): os.remove(out)
subs = sorted(glob.glob(f"{sd}/s*.tre"))
env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(TM_PY) + "/../lib")
t0 = time.time()
p = subprocess.run([TM_PY, f"{TM_DIR}/treemerge.py", "-s", f"{d}/guide.tre", "-m", mat, "-x", mat + "_taxlist",
                    "-o", out, "-p", PAUP, "-w", work, "-t", *subs], capture_output=True, text=True, env=env)
sec = time.time() - t0
if p.returncode != 0 or not os.path.exists(out):
    print(p.stdout[-2000:], p.stderr[-3000:]); sys.exit(1)
T = read_tree(f"{d}/true.tre")
from phylo import is_induced
M = read_tree(out)
res = dict(dir=d, subtrees=st, fn_treemerge=fn_fp(M, T)[0], ok=all(is_induced(M, read_tree(s)) for s in subs), sec=sec)
json.dump(res, open(f"{sd}/treemerge.json", "w")); print(json.dumps(res), flush=True)
