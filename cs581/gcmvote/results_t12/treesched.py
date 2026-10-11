# AI-assisted (Claude), exploration code for CS581 project
# Starts FastTree (protbench trees.py, one method per process) as alignments finish; max 4 at once.
import filecmp, json, os, subprocess, time
R = "/opt/work/gcmtrees/reps/SIMHIGH_R12"; V = "/opt/work/gcmvote/reps/SIMHIGH_R12"
T = "/opt/work/gcmtrees/trees/SIMHIGH_R12_d0"; OUT = "/opt/work/gcmtrees/trees_R12.jsonl"
PY = "/opt/mm/root/envs/pasta183/bin/python"; TP = "/home/user/scratch/cs581/protbench/code/trees.py"
TRUE = "/opt/data/sim/SIMHIGH/R12/tree.nwk"
gg = {"magus": "linsi", "es4": "linsi#es4", "recipe": "wsoft0.03:linsi&fftns2#es4", "es3": "linsi#es3",
      "hard": "linsi&fftns2-op3"}
ggdir = {"linsi": "linsi", "linsi#es4": "linsi_es_es4", "wsoft0.03:linsi&fftns2#es4": "wsoft0.03_c_linsi_i_fftns2_es_es4",
         "linsi#es3": "linsi_es_es3", "linsi&fftns2-op3": "linsi_i_fftns2-op3"}
vote = ["hard", "soft", "soft2", "soft4", "hard-bb", "soft-bb", "hard+mask"]
def rows(p, key):
    return {key(json.loads(l)) for l in open(p)} if os.path.exists(p) else set()
def ready():
    out = {"true": os.path.join(T, "true.fasta")}
    g = rows(R + "/results.jsonl", lambda r: r["variant"])
    for m, v in gg.items():
        if v in g: out[m] = "{}/variants/{}/out.fasta".format(R, ggdir[v])
    vv = rows(V + "/results.jsonl", lambda r: (r["variant"], r["B"]))
    for v in vote:
        if (v, 10) in vv:
            out["vote_" + v] = "{}/vote/{}/out{}.fasta".format(V, v, ".masked" if v.endswith("+mask") else "")
    return out
order = ["true", "magus", "es4", "recipe", "es3", "hard"] + ["vote_" + v for v in vote]
started, procs, reuse = {}, {}, {}
done = rows(OUT, lambda r: r["method"])
for m in done: started[m] = None
started["true"] = os.path.join(T, "true.fasta")
while True:
    for m, p in list(procs.items()):
        if p.poll() is not None: del procs[m]; print("finished", m, p.returncode, flush=True)
    rd = ready()
    for m in order:
        if m in started or m in reuse or m not in rd or len(procs) >= 4: continue
        src = rd[m]
        twin = next((o for o, s in started.items() if s and filecmp.cmp(src, s, shallow=False)), None)
        if twin:
            reuse[m] = twin; print("reuse", m, "=", twin, flush=True); continue
        dst = os.path.join(T, m + ".fasta")
        if os.path.lexists(dst): os.remove(dst)
        os.symlink(src, dst); started[m] = src
        procs[m] = subprocess.Popen(["nice", "-n", "5", PY, TP, T, TRUE, OUT, m],
                                    stdout=open("/tmp/claude-0/logs/tree_{}.log".format(m), "w"), stderr=subprocess.STDOUT)
        print("start", m, flush=True)
    json.dump(reuse, open("/opt/work/gcmtrees/trees_R12_reuse.json", "w"))
    if not procs and all(m in started or m in reuse for m in order): break
    time.sleep(30)
print("TREES_DONE", flush=True)
