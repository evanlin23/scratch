"""MAGUS-fast pilot driver (restartable; one step at a time on an idle 4-core machine).

    python mf_bench.py STEP NAME [NAME ...]

STEP:
  default   MAGUS end to end, paper flags, pure MAGUS (Python graph builder, 1-thread MCL), profiled
            (mfrun.py MF_PROF); keeps subalignments + backbones under WORK/NAME/inputs
  graph     merge-only on those inputs, 10 x 200 backbones: Python builder vs vectorized builder vs
            vectorized + MCL -te 4; graphs, clusterings and outputs compared for identity
  bb100     size-100 backbones: per default backbone set, a seeded half per subset (4 of 8), realigned
            with MAGUS's own L-INS-i command, 4 at a time with --thread 4 (as MAGUS schedules them)
  sweep     merge-only N in {10,6,4,3} x s in {200,100} x pruning {off, K=ceil(0.4N)}, vectorized graph +
            MCL -te 4; first N backbones of each size
  e2e:SET   end-to-end MAGUS with setting SET (e.g. n4s100p1), vectorized graph, MCL -te 4, profiled
  tree      FastTree on the true alignment and on every scored alignment kept for NAME (simulated sets)

State: WORK/NAME/state.json (atomic writes); a finished step is skipped on rerun.
"""
import concurrent.futures
import json
import math
import os
import random
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "cs581", "code"))
from gcmx import fasta, score  # noqa: E402
from gcmx.e2e_bench import magus_flags, seq_type  # noqa: E402

WORK = os.environ.get("MF_WORK", "/opt/work/mf")
RES = os.path.join(REPO, "cs581", "magusfast", "results")
T = "4"
MERGE_FLAGS = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
               "-f", "4"]
DATA = {
    "BBA0101_R0": ("cs581/data/balibase_clean/RV100_BBA0101.fasta", None),
    "BBA0190_R0": ("cs581/data/balibase_clean/RV100_BBA0190.fasta", None),
    "BBA0067_R0": ("cs581/data/balibase_clean/RV100_BBA0067.fasta", None),
    "1000M2_R0": ("/opt/data/Datasets/ROSE/1000M2/R0/rose.aln.true.fasta", "/opt/data/Datasets/ROSE/1000M2/R0/rose.tt"),
    "1000L1_R0": ("/opt/data/Datasets/ROSE/1000L1/R0/rose.aln.true.fasta", "/opt/data/Datasets/ROSE/1000L1/R0/rose.tt"),
    "1000M3_R0": ("/opt/data/Datasets/ROSE/1000M3/R0/rose.aln.true.fasta", "/opt/data/Datasets/ROSE/1000M3/R0/rose.tt"),
    "RNASim1000_R0": ("/opt/data/Datasets/RNASim/1000/R0/true_align.txt", None),
    "RNASim1000_R1": ("/opt/data/Datasets/RNASim/1000/R1/true_align.txt", None),
    "SIMHIGH_R1": ("/opt/data/sim/SIMHIGH/R1/sim.fa", "/opt/data/sim/SIMHIGH/R1/tree.nwk"),
    "SIMHIGH_R2": ("/opt/data/sim/SIMHIGH/R2/sim.fa", "/opt/data/sim/SIMHIGH/R2/tree.nwk"),
}


def K_of(n):
    return math.ceil(0.4 * n)


def setting_name(n, s, p):
    return "n{}s{}p{}".format(n, s, p)


def parse_setting(name):
    n, rest = name[1:].split("s")
    s, p = rest.split("p")
    return int(n), int(s), int(p)


class Run:
    def __init__(self, name):
        self.name = name
        src, self.tree = DATA[name]
        src = src if os.path.isabs(src) else os.path.join(REPO, src)
        self.w = os.path.join(WORK, name)
        os.makedirs(self.w, exist_ok=True)
        self.true = os.path.join(self.w, "true.fasta")
        self.unaligned = os.path.join(self.w, "unaligned.fasta")
        if not os.path.exists(self.unaligned):
            ref = fasta.upper(fasta.read(src))
            fasta.write(ref, self.true)
            fasta.write(fasta.ungap(ref), self.unaligned)
        self.state_path = os.path.join(self.w, "state.json")
        self.state = json.load(open(self.state_path)) if os.path.exists(self.state_path) else {
            "dataset": name, "datatype": seq_type(fasta.read(self.true)), "nseq": len(fasta.read(self.true))}
        self.inputs = os.path.join(self.w, "inputs")

    def save(self, key, data):
        self.state[key] = data
        with open(self.state_path + ".tmp", "w") as f:
            json.dump(self.state, f)
        os.replace(self.state_path + ".tmp", self.state_path)
        os.makedirs(os.path.join(RES, "state"), exist_ok=True)
        shutil.copy(self.state_path, os.path.join(RES, "state", self.name + ".json"))
        print(json.dumps({self.name: {key: data}})[:600], flush=True)

    def acc(self, path):
        s = score.fastsp(self.true, path)
        return {"SPFN": s["SPFN"], "SPFP": s["SPFP"], "err": round(50 * (s["SPFN"] + s["SPFP"]), 3)}


def timed(cmd, log, env=None):
    start = time.time()
    with open(log, "w") as f:
        p = subprocess.Popen(cmd, cwd=HERE, stdout=f, stderr=subprocess.STDOUT, env=env)
        _, status, ru = os.wait4(p.pid, 0)
    if os.waitstatus_to_exitcode(status) != 0:
        raise RuntimeError("failed: " + " ".join(cmd) + " (log " + log + ")")
    return round(time.time() - start, 1), round(ru.ru_utime + ru.ru_stime, 1)


def mfrun(args, log, prof=None, esk=1):
    env = dict(os.environ)
    if prof:
        env["MF_PROF"] = prof
    env["MF_ESK"] = str(esk)
    return timed([sys.executable, os.path.join(HERE, "mfrun.py")] + args, log, env)


def summarize_prof(prof):
    """Stage split of one profiled run (seconds since start)."""
    p = json.load(open(prof))
    st = p["stages"]
    tasks = p["tasks"]
    span = lambda k: sum(e - s for s, e in st.get(k, []))  # noqa: E731
    pool = [t for t in tasks if t["kind"] in ("subset", "backbone")]
    out = {"total": p["total"]}
    dec_end = st["decomposition"][0][1] if "decomposition" in st else 0.0
    out["decomposition"] = round(dec_end, 1)
    if pool:
        pool_end = max(t["end"] for t in pool)
        out["pool"] = round(pool_end - min(t["start"] for t in pool), 1)
        out["pool_end"] = round(pool_end, 1)
        for kind in ("subset", "backbone"):
            ts = [t for t in tasks if t["kind"] == kind]
            out[kind + "_cpu"] = round(sum(t["cpu"] for t in ts), 1)
            out[kind + "_n"] = len(ts)
            out[kind + "_end"] = round(max((t["end"] for t in ts), default=0), 1)
            out[kind + "_task_wall"] = round(sum(t["end"] - t["start"] for t in ts), 1)
        out["bb_cpu_each"] = [t["cpu"] for t in sorted((t for t in tasks if t["kind"] == "backbone"),
                                                        key=lambda t: int(os.path.basename(t["out"]).split("_")[1]))]
    else:
        pool_end = 0.0
    bm = st.get("buildMatrix", [[0, 0]])[0]
    out["graph_compute_python"] = round(span("graph_add_python"), 1)
    out["graph_exposed"] = round(bm[1] - max(bm[0], pool_end), 1)  # build time after the last MAFFT task
    out["graph_buildMatrix"] = round(bm[1] - bm[0], 1)
    out["graph_write"] = round(span("writeGraph"), 1)
    out["mcl"] = round(span("cluster"), 1)
    out["trace"] = round(span("trace") + span("optimize"), 1)
    out["write"] = round(span("write"), 1)
    out["mcl_cpu"] = round(sum(t["cpu"] for t in tasks if t["kind"] == "mcl"), 1)
    return out


def step_default(r):
    if "default" in r.state and os.path.isdir(r.inputs):
        return
    d = os.path.join(r.w, "magus")
    shutil.rmtree(d, ignore_errors=True)
    shutil.rmtree(r.inputs, ignore_errors=True)
    out = os.path.join(r.w, "default.fasta")
    if os.path.exists(out):
        os.remove(out)
    prof = os.path.join(r.w, "default.prof.json")
    wall, cpu = mfrun(["--gcmx-fastgraph", "false", "-np", T, "-d", d, "-i", r.unaligned, "-o", out] + magus_flags(25),
                      os.path.join(r.w, "default.log"), prof)
    shutil.copytree(os.path.join(d, "subalignments"), os.path.join(r.inputs, "subalignments"))
    os.makedirs(os.path.join(r.inputs, "backbones"))
    for f in os.listdir(os.path.join(d, "graph")):
        if f.startswith("backbone_") and f.endswith(("_unalign.txt", "_mafft.txt")):
            shutil.copy(os.path.join(d, "graph", f), os.path.join(r.inputs, "backbones"))
    shutil.copy(prof, os.path.join(RES, "prof", r.name + ".default.json")) if os.path.isdir(os.path.join(RES, "prof")) else None
    shutil.rmtree(d, ignore_errors=True)
    r.save("default", {"wall": wall, "cpu": cpu, **r.acc(out), "prof": summarize_prof(prof)})


def bb_dir(r, n, s):
    """Directory holding exactly the first n backbones of size s."""
    d = os.path.join(r.w, "bbsel", "n{}s{}".format(n, s))
    if os.path.isdir(d):
        return d
    src = os.path.join(r.inputs, "backbones" if s == 200 else "bb100")
    os.makedirs(d + ".tmp", exist_ok=True)
    for b in range(1, n + 1):
        shutil.copy(os.path.join(src, "backbone_{}_mafft.txt".format(b)), d + ".tmp")
    os.replace(d + ".tmp", d)
    return d


def merge(r, label, n, s, esk, fastgraph=True, mclthreads=4, keep_graph=False):
    d = os.path.join(r.w, "m_" + label)
    shutil.rmtree(d, ignore_errors=True)
    out = os.path.join(r.w, label + ".fasta")
    if os.path.exists(out):
        os.remove(out)
    prof = os.path.join(r.w, label + ".prof.json")
    args = ["--gcmx-fastgraph", "true" if fastgraph else "false", "-np", T, "-d", d,
            "-s", os.path.join(r.inputs, "subalignments"), "-b", bb_dir(r, n, s), "-o", out] + MERGE_FLAGS
    if mclthreads > 1:
        args = ["--gcmx-mclthreads", str(mclthreads)] + args
    wall, cpu = mfrun(args, os.path.join(r.w, label + ".log"), prof, esk)
    res = {"merge_wall": wall, "merge_cpu": cpu, **r.acc(out), "prof": summarize_prof(prof)}
    if not keep_graph:
        shutil.rmtree(d, ignore_errors=True)
    return res, d, out


def graph_entries(path):
    import numpy as np
    a = np.loadtxt(path, dtype=np.int64)
    a = a[np.lexsort((a[:, 1], a[:, 0]))]
    return a


def step_graph(r):
    if "graph" in r.state:
        return
    py, dpy, opy = merge(r, "g_python", 10, 200, 1, fastgraph=False, mclthreads=1, keep_graph=True)
    fg, dfg, ofg = merge(r, "g_fast", 10, 200, 1, fastgraph=True, mclthreads=1, keep_graph=True)
    fm, dfm, ofm = merge(r, "g_fast_mcl4", 10, 200, 1, fastgraph=True, mclthreads=4, keep_graph=True)
    import numpy as np
    gp, gf = graph_entries(os.path.join(dpy, "graph", "graph.txt")), graph_entries(os.path.join(dfg, "graph", "graph.txt"))
    def clusters(x):
        return sorted(tuple(sorted(map(int, l.split()))) for l in open(os.path.join(x, "graph", "clusters.txt")) if l.strip())

    def same_homologies(a, b):  # identical aligned residue pairs and length (column order may differ)
        s = score.fastsp(a, b)
        return s["SPFN"] == 0 and s["SPFP"] == 0 and s["LenRef"] == s["LenEst"]
    cl = [clusters(x) for x in (dpy, dfg, dfm)]
    ident = {"graph_identical": bool(gp.shape == gf.shape and (gp == gf).all()), "graph_entries": int(len(gp)),
             "clusters_identical_fast": cl[0] == cl[1], "clusters_identical_mcl4": cl[0] == cl[2],
             "output_identical_fast": same_homologies(opy, ofg), "output_identical_mcl4": same_homologies(opy, ofm),
             "python_merge_reproduces_default": same_homologies(os.path.join(r.w, "default.fasta"), opy)}
    for x in (dpy, dfg, dfm):
        shutil.rmtree(x, ignore_errors=True)
    r.save("graph", {"python": py, "fast": fg, "fast_mcl4": fm, **ident})


def subset_of(r):
    m = {}
    sd = os.path.join(r.inputs, "subalignments")
    for f in sorted(os.listdir(sd)):
        for t in fasta.read(os.path.join(sd, f)):
            m[t] = f
    return m


def mafft_bin():
    from magus.configuration import Configs
    return Configs.mafftPath


def step_bb100(r, nbb=10):
    if "bb100" in r.state:
        return
    dst = os.path.join(r.inputs, "bb100")
    os.makedirs(dst, exist_ok=True)
    sub = subset_of(r)
    unal = fasta.read(r.unaligned)
    todo = []
    for b in range(1, nbb + 1):
        seqs = fasta.read(os.path.join(r.inputs, "backbones", "backbone_{}_unalign.txt".format(b)))
        groups = {}
        for t in seqs:
            groups.setdefault(sub[t], []).append(t)
        rng = random.Random(1000 + b)
        keep = []
        for g in sorted(groups):
            ts = sorted(groups[g])
            rng.shuffle(ts)
            keep += ts[:max(1, len(ts) // 2)]
        u = os.path.join(dst, "backbone_{}_unalign.txt".format(b))
        fasta.write({t: unal[t] for t in keep}, u)
        todo.append(b)
    times = {}

    def one(b):
        u = os.path.join(dst, "backbone_{}_unalign.txt".format(b))
        out = os.path.join(dst, "backbone_{}_mafft.txt".format(b))
        if os.path.exists(out):
            return b, None
        cmd = [mafft_bin(), "--localpair", "--maxiterate", "1000", "--ep", "0.123", "--quiet", "--thread", T,
               "--anysymbol", u]
        s = time.time()
        with open(out + ".tmp", "w") as o:
            p = subprocess.Popen(cmd, stdout=o, stderr=subprocess.DEVNULL)
            _, status, ru = os.wait4(p.pid, 0)
        assert os.waitstatus_to_exitcode(status) == 0
        os.replace(out + ".tmp", out)
        return b, {"wall": round(time.time() - s, 1), "cpu": round(ru.ru_utime + ru.ru_stime, 1), "nseq": len(keep)}

    start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(T)) as pool:
        for b, t in pool.map(one, todo):
            times[b] = t
    r.save("bb100", {"wall": round(time.time() - start, 1),
                     "cpu_each": [times[b]["cpu"] if times[b] else None for b in todo],
                     "nseq_each": [len(fasta.read(os.path.join(dst, "backbone_{}_unalign.txt".format(b)))) for b in todo]})


def step_sweep(r):
    sw = r.state.get("sweep", {})
    for s in (200, 100):
        for n in (10, 6, 4, 3):
            for p in (0, 1):
                name = setting_name(n, s, p)
                if name in sw:
                    continue
                res, _, out = merge(r, "sw_" + name, n, s, K_of(n) if p else 1)
                os.replace(out, os.path.join(r.w, "sw_" + name + ".fasta"))
                sw[name] = res
                r.save("sweep", sw)


def step_e2e(r, setting):
    key = "e2e_" + setting
    if key in r.state:
        return
    n, s, p = parse_setting(setting)
    d = os.path.join(r.w, key)
    shutil.rmtree(d, ignore_errors=True)
    out = os.path.join(r.w, key + ".fasta")
    if os.path.exists(out):
        os.remove(out)
    prof = os.path.join(r.w, key + ".prof.json")
    flags = magus_flags(25)
    flags[flags.index("-r") + 1] = str(n)
    flags[flags.index("-m") + 1] = str(s)
    wall, cpu = mfrun(["--gcmx-fastgraph", "true", "--gcmx-mclthreads", T, "-np", T, "-d", d, "-i", r.unaligned,
                       "-o", out] + flags, os.path.join(r.w, key + ".log"), prof, K_of(n) if p else 1)
    shutil.rmtree(d, ignore_errors=True)
    r.save(key, {"wall": wall, "cpu": cpu, **r.acc(out), "prof": summarize_prof(prof)})


def nrf(true_tree, est_tree):
    import dendropy
    from dendropy.calculate import treecompare
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(path=true_tree, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    t2 = dendropy.Tree.get(path=est_tree, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    t1.is_rooted = t2.is_rooted = False
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    fp, fn = treecompare.false_positives_and_negatives(t1, t2)
    nint = len(t1.leaf_nodes()) - 3
    return round(100 * (fp + fn) / (2 * nint), 2)


def step_tree(r):
    if not r.tree:
        return
    tr = r.state.get("tree", {})
    prot = r.state["datatype"] == "protein"
    alns = {"true": r.true, "default": os.path.join(r.w, "default.fasta")}
    for f in sorted(os.listdir(r.w)):
        if f.endswith(".fasta") and (f.startswith("e2e_") or f.startswith("sw_")):
            alns[f[:-6]] = os.path.join(r.w, f)
    only = os.environ.get("MF_TREES")  # comma list of labels to restrict to
    for label, path in alns.items():
        if label in tr or (only and label not in only.split(",") and label not in ("true", "default")):
            continue
        out = os.path.join(r.w, "tree_" + label + ".nwk")
        cmd = ["FastTree", "-quiet", "-nopr"] + (["-lg", "-gamma"] if prot else ["-nt", "-gtr", "-gamma"]) + [path]
        with open(out, "w") as o:
            subprocess.run(cmd, stdout=o, stderr=subprocess.DEVNULL, check=True)
        tr[label] = nrf(r.tree, out)
        r.save("tree", tr)


def main():
    step, names = sys.argv[1], sys.argv[2:]
    os.makedirs(os.path.join(RES, "prof"), exist_ok=True)
    for name in names:
        r = Run(name)
        if step == "default":
            step_default(r)
        elif step == "graph":
            step_graph(r)
        elif step == "bb100":
            step_bb100(r)
        elif step == "sweep":
            step_sweep(r)
        elif step.startswith("e2e:"):
            step_e2e(r, step[4:])
        elif step == "tree":
            step_tree(r)
        else:
            raise SystemExit("unknown step " + step)


if __name__ == "__main__":
    main()
