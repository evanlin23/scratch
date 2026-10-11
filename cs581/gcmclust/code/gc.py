"""Swap the clustering step of MAGUS's GCM (merge-only, paired vs MAGUS on identical subsets and backbones).

    python3 gc.py prep  NAME REP_DIR TRUE INPUTS     # INPUTS = inputs.tar.xz or a dir holding inputs/
    python3 gc.py build REP_DIR                      # MAGUS control merge + graph.npz (weights, backbone support)
    python3 gc.py run   REP_DIR VARIANT [...]        # -> REP_DIR/results.jsonl (restartable, lock per variant)

VARIANT = GRAPH:METHOD[:PARAM]
  GRAPH   raw  MAGUS's graph;  es4  cross-subset edges supported by < 4 of the 10 backbones deleted
          (gcmgen's `linsi#es4`; the trace step also sees the filtered graph, as in gcmgen)
  METHOD  mcl:I       MCL (MAGUS's bundled binary, same abc input as MAGUS) with inflation I; raw:mcl:4 = MAGUS
          leidmod:g   Leiden, modularity (RBConfiguration) with resolution g
          leidcpm:g   Leiden, CPM with resolution g on degree-normalised weights w/sqrt(s_a s_b)
          louvain:g   Louvain (igraph multilevel), modularity with resolution g
          cc:t        connected components of the cross-subset edges supported by >= t backbones
          agglo:t     greedy constrained agglomeration: cross-subset edges with support >= t in decreasing
                      order of normalised weight; two clusters are joined only if they share no subset
          lpa         label propagation (igraph)
Only the clustering changes; MAGUS's own purge + minclusters trace runs on the clusters (graph/clusters.txt
is pre-written, which MAGUS reads instead of running MCL).
"""
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, "..", "..", "code"))
sys.path.insert(0, CODE)
from gcmx import fasta, score  # noqa: E402

FLAGS = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
         "-f", "4"]
NP = os.environ.get("GC_NP", "1")
TIMEOUT = 1800  # merge wall limit (s); a variant that exceeds it on any replicate is recorded as failed


def mcl_bin():
    return os.path.join(CODE, "MAGUS", "magus", "tools", "mcl", "mcl")


# ------------------------------------------------------------------ replicate preparation

def prep(name, rep, true, inputs):
    if os.path.exists(os.path.join(rep, "true.fasta")):
        return
    os.makedirs(rep, exist_ok=True)
    if inputs.endswith(".tar.xz"):
        subprocess.run(["tar", "xJf", inputs, "-C", rep], check=True)
    else:
        shutil.copytree(os.path.join(inputs, "inputs"), os.path.join(rep, "inputs"))
    for f in os.listdir(os.path.join(rep, "inputs", "backbones")):
        if "unalign" in f:  # fresh bbtool_bench draws keep MAGUS's unaligned backbone files next to the aligned ones
            os.remove(os.path.join(rep, "inputs", "backbones", f))
    for d in ("subalignments", "backbones"):
        p = os.path.join(rep, "inputs", d)
        for f in os.listdir(p):
            fasta.write(fasta.upper(fasta.read(os.path.join(p, f))), os.path.join(p, f))
    fasta.write(fasta.upper(fasta.read(true)), os.path.join(rep, "true.fasta"))
    json.dump({"name": name, "true": true, "inputs": inputs}, open(os.path.join(rep, "prep.json"), "w"))


def acc(rep, out):
    s = score.fastsp(os.path.join(rep, "true.fasta"), out)
    return {"SPFN": s["SPFN"], "SPFP": s["SPFP"], "err": round(100 * (s["SPFN"] + s["SPFP"]) / 2, 4),
            "LenEst": s["LenEst"], "LenRef": s["LenRef"]}


def magus(rep, work, out, extra_pre=(), script=None):
    argv = [sys.executable] + ([script] + list(extra_pre) if script else ["-m", "gcmx.run_magus"])
    argv += ["--gcmx-fastgraph", "false", "-np", NP, "-d", work, "-s", os.path.join(rep, "inputs", "subalignments"),
             "-b", os.path.join(rep, "inputs", "backbones"), "-o", out] + FLAGS
    start = time.time()
    with open(os.path.join(os.path.dirname(out), "magus.log"), "w") as log:
        subprocess.run(argv, cwd=CODE, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=TIMEOUT)
    return round(time.time() - start, 1)


def build(rep):
    """Control merge (MAGUS itself) and graph.npz; keeps graph_raw.txt (MAGUS's abc graph) for MCL."""
    if os.path.exists(os.path.join(rep, "graph.npz")):
        return
    d = os.path.join(rep, "ctrl")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    out = os.path.join(d, "out.fasta")
    wall = magus(rep, os.path.join(d, "work"), out, [os.path.join(rep, "graph.tmp.npz")],
                 os.path.join(HERE, "build_graph.py"))
    shutil.move(os.path.join(d, "work", "graph", "graph.txt"), os.path.join(rep, "graph_raw.txt"))
    shutil.rmtree(os.path.join(d, "work"))
    row = {"rep": os.path.basename(rep), "variant": "magus", "merge_wall": wall, **acc(rep, out)}
    json.dump(row, open(os.path.join(rep, "ctrl.json"), "w"))
    os.replace(os.path.join(rep, "graph.tmp.npz"), os.path.join(rep, "graph.npz"))
    print(json.dumps(row), flush=True)


# ------------------------------------------------------------------ graphs and clusterers

class Graph:
    def __init__(self, rep, kind):
        z = np.load(os.path.join(rep, "graph.npz"))
        a, b, w, nbb, self.sub = z["a"], z["b"], z["w"].astype(np.float64), z["nbb"], z["sub"]
        self.n = len(self.sub)
        if kind == "es4":
            keep = (a == b) | (nbb >= 4)  # exactly gcmgen/protcons run_edgesup K=4
        elif kind == "raw":
            keep = np.ones(len(a), bool)
        else:
            raise ValueError(kind)
        self.a, self.b, self.w, self.nbb = a[keep], b[keep], w[keep], nbb[keep]
        self.kind = kind

    def write_abc(self, path):
        """The graph as MAGUS writes it (directed abc lines, both directions, self loops)."""
        with open(path, "w") as f:
            for x, y, w in zip(self.a.tolist(), self.b.tolist(), self.w.astype(np.int64).tolist()):
                f.write("{} {} {}\n".format(x, y, w))

    def undirected(self, cross_only=False, min_nbb=1):
        """a < b edges (no self loops), optionally cross-subset only and with support >= min_nbb."""
        sel = self.a < self.b
        if cross_only:
            sel &= self.sub[self.a] != self.sub[self.b]
        if min_nbb > 1:
            sel &= self.nbb >= min_nbb
        return self.a[sel], self.b[sel], self.w[sel], self.nbb[sel]

    def strength(self):
        a, b, w, _ = self.undirected()
        s = np.bincount(a, w, self.n) + np.bincount(b, w, self.n)
        return np.maximum(s, 1e-9)


def to_clusters(labels):
    order = np.argsort(labels, kind="stable")
    lab = labels[order]
    cuts = np.flatnonzero(np.diff(lab)) + 1
    return [g.tolist() for g in np.split(order, cuts)]


def igraph_of(G):
    import igraph as ig
    a, b, w, _ = G.undirected()
    g = ig.Graph(n=G.n, edges=np.column_stack([a, b]).tolist(), directed=False)
    return g, w


def cluster(G, method, param, work):
    if method == "mcl":
        abc = os.path.join(work, "graph.txt")
        out = os.path.join(work, "mcl.out")
        subprocess.run([mcl_bin(), abc, "--abc", "-o", out, "-I", str(param), "-te", "1"], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return [[int(t) for t in line.split()] for line in open(out) if line.strip()]
    if method in ("leidmod", "leidcpm"):
        import leidenalg as la
        g, w = igraph_of(G)
        if method == "leidmod":
            p = la.find_partition(g, la.RBConfigurationVertexPartition, weights=w.tolist(),
                                  resolution_parameter=float(param), seed=1)
        else:
            s = G.strength()
            a, b, _, _ = G.undirected()
            wn = w / np.sqrt(s[a] * s[b])
            p = la.find_partition(g, la.CPMVertexPartition, weights=wn.tolist(),
                                  resolution_parameter=float(param), seed=1)
        return [list(c) for c in p]
    if method == "louvain":
        import random
        random.seed(1)
        g, w = igraph_of(G)
        p = g.community_multilevel(weights=w.tolist(), resolution=float(param))
        return [list(c) for c in p]
    if method == "lpa":
        import random
        random.seed(1)
        g, w = igraph_of(G)
        return [list(c) for c in g.community_label_propagation(weights=w.tolist())]
    if method == "cc":
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        a, b, w, _ = G.undirected(cross_only=True, min_nbb=int(param))
        m = coo_matrix((np.ones(len(a)), (a, b)), shape=(G.n, G.n))
        _, lab = connected_components(m, directed=False)
        return to_clusters(lab)
    if method == "agglo":
        s = G.strength()
        a, b, w, _ = G.undirected(cross_only=True, min_nbb=int(param))
        order = np.argsort(-(w / np.sqrt(s[a] * s[b])), kind="stable")
        parent = list(range(G.n))
        subs = [{int(x)} for x in G.sub]

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for i in order.tolist():
            x, y = find(int(a[i])), find(int(b[i]))
            if x == y or not subs[x].isdisjoint(subs[y]):
                continue
            if len(subs[x]) < len(subs[y]):
                x, y = y, x
            parent[y] = x
            subs[x] |= subs[y]
            subs[y] = None
        return to_clusters(np.array([find(x) for x in range(G.n)]))
    raise ValueError(method)


def stats(G, clusters):
    sizes = np.array([len(c) for c in clusters])
    multi = sizes[sizes > 1]
    viol = sum(1 for c in clusters if len(c) > 1 and len(set(G.sub[c].tolist())) < len(c))
    return {"n_clusters": int(len(multi)), "mean_size": round(float(multi.mean()), 2) if len(multi) else 0,
            "max_size": int(sizes.max()), "frac_singleton": round(float((sizes == 1).sum() / G.n), 4),
            "frac_violating": round(viol / max(len(multi), 1), 4), "nodes": int(G.n)}


def run(rep, names):
    res = os.path.join(rep, "results.jsonl")
    done = {json.loads(l)["variant"] for l in open(res)} if os.path.exists(res) else set()
    for name in names:
        if name in done:
            continue
        vd = os.path.join(rep, "variants", name.replace(":", "_"))
        lock = vd + ".lock"
        os.makedirs(os.path.dirname(lock), exist_ok=True)
        try:
            os.close(os.open(lock, os.O_CREAT | os.O_EXCL))
        except FileExistsError:
            continue
        if os.path.exists(res) and name in {json.loads(l)["variant"] for l in open(res)}:
            os.remove(lock)  # another lane finished it meanwhile
            continue
        try:
            kind, method, *param = name.split(":")
            param = param[0] if param else None
            shutil.rmtree(vd, ignore_errors=True)
            gdir = os.path.join(vd, "work", "graph")
            os.makedirs(gdir)
            t0 = time.time()
            G = Graph(rep, kind)
            if kind == "raw":
                shutil.copy(os.path.join(rep, "graph_raw.txt"), os.path.join(gdir, "graph.txt"))
            else:
                G.write_abc(os.path.join(gdir, "graph.txt"))
            t1 = time.time()
            cl = cluster(G, method, param, gdir)
            t2 = time.time()
            st = stats(G, cl)
            with open(os.path.join(gdir, "clusters.txt"), "w") as f:
                for c in cl:
                    if len(c) > 1:
                        f.write(" ".join(map(str, c)) + "\n")
            out = os.path.join(vd, "out.fasta")
            m_wall = magus(rep, os.path.join(vd, "work"), out)
            ntrace = sum(1 for _ in open(os.path.join(gdir, "trace.txt")))
            row = {"rep": os.path.basename(rep), "variant": name, "graph_s": round(t1 - t0, 1),
                   "cluster_s": round(t2 - t1, 1), "merge_wall": m_wall, **st, "trace_clusters": ntrace,
                   **acc(rep, out)}
            shutil.rmtree(os.path.join(vd, "work"), ignore_errors=True)
            with open(res, "a") as f:
                f.write(json.dumps(row) + "\n")
            print(json.dumps(row), flush=True)
        except Exception as e:  # keep the lane going; a failure row is recorded (timeouts count as failures)
            print("FAILED", rep, name, repr(e), flush=True)
            with open(res, "a") as f:
                f.write(json.dumps({"rep": os.path.basename(rep), "variant": name, "failed": repr(e)[:200]}) + "\n")
            shutil.rmtree(os.path.join(vd, "work"), ignore_errors=True)
            os.remove(lock)
        else:
            os.remove(lock)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "prep":
        prep(sys.argv[2], os.path.abspath(sys.argv[3]), sys.argv[4], sys.argv[5])
    elif cmd == "build":
        build(os.path.abspath(sys.argv[2]))
    elif cmd == "run":
        run(os.path.abspath(sys.argv[2]), sys.argv[3:])
