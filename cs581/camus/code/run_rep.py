"""Run CAMUS variants on one replicate and score them against the true network.

usage: python3 run_rep.py COND REP GT OUTDIR [--variants a,b,...]
  COND  n15|n25|n50|...   REP 00..49   GT g_500|iqtree_500|g_true
Writes OUTDIR/COND_GT_REP.jsonl (one line per variant; restartable: done variants skipped).
"""
import csv
import json
import os
import subprocess
import sys
import time
import warnings

import treeswift as ts

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import netmetric as nm  # noqa: E402
import quartets as qq  # noqa: E402

DATA = "/opt/data/camus/sim/camus-dataset"
PUB = "/opt/data/camus/inf/inferred-networks"
BIN = "/opt/mm/root/envs/phy/bin"
CAMUS = "/opt/tools/camus/camus-z"
NPROC = int(os.environ.get("CAMUS_NPROC", "1"))

# name -> (base tree, t, z, zmode)
VARIANTS = {
    "default": ("astral_pub", 0.5, 0, ""),
    "t0": ("astral_pub", 0.0, 0, ""),
    "t0.2": ("astral_pub", 0.2, 0, ""),
    "t0.3": ("astral_pub", 0.3, 0, ""),
    "t0.7": ("astral_pub", 0.7, 0, ""),
    "t0.8": ("astral_pub", 0.8, 0, ""),
    "z2and": ("astral_pub", 0.5, 2.0, "and"),
    "z3and": ("astral_pub", 0.5, 3.0, "and"),
    "z3only": ("astral_pub", 0.0, 3.0, "only"),
    "z5only": ("astral_pub", 0.0, 5.0, "only"),
    "z3or": ("astral_pub", 0.5, 3.0, "or"),
    "tqmc": ("tqmc", 0.5, 0, ""),
    "wastral": ("wastral", 0.5, 0, ""),
    "true_major": ("true_major", 0.5, 0, ""),
    "swap": ("swap", 0.5, 0, ""),  # other displayed tree of default's 1-reticulation network
}


def root_out(newick):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = ts.read_tree_newick(newick)
        if "OUT" not in [c.label for c in t.root.children]:
            out = [x for x in t.traverse_leaves() if x.label == "OUT"][0]
            out.edge_length = 1.0
            t.reroot(out, length=0.5)
        t.suppress_unifurcations()
        t.resolve_polytomies()
        for n in t.traverse_preorder():
            n.edge_length = None
            if not n.is_leaf():
                n.label = None
        t.is_rooted = True
    return t.newick().replace("[&R] ", "")


def run_camus(base_nwk, genes, t, z, zmode, prefix):
    with open(prefix + ".base.nwk", "w") as f:
        f.write(base_nwk + "\n")
    env = dict(os.environ, CAMUS_Z=str(z), CAMUS_ZMODE=zmode)
    t0 = time.time()
    p = subprocess.run([CAMUS, "-n", str(NPROC), "-t", str(t), "-o", prefix, prefix + ".base.nwk", genes],
                       env=env, capture_output=True, text=True)
    dt = time.time() - t0
    if p.returncode != 0:
        raise RuntimeError(p.stderr[-2000:])
    rows = list(csv.reader(open(prefix + ".csv")))[1:]
    nets = {int(r[0]): (r[2], float(r[1])) for r in rows}
    nq = None
    for line in open(prefix + ".log"):
        if "containing" in line and "quartets not in the constraint tree" in line:
            nq = int(line.split("containing")[1].split()[0])
    return nets, dt, nq


def main():
    cond, rep, gt, outdir = sys.argv[1:5]
    want = list(VARIANTS)
    if "--variants" in sys.argv:
        want = sys.argv[sys.argv.index("--variants") + 1].split(",")
    os.makedirs(outdir, exist_ok=True)
    work = f"/opt/runs/camus/{cond}_{gt}_{rep}"
    os.makedirs(work, exist_ok=True)
    out = os.path.join(outdir, f"{cond}_{gt}_{rep}.jsonl")
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)["variant"] for l in open(out)}
    rd = f"{DATA}/{cond}/{rep}"
    genes = f"{rd}/{gt}.nwk"
    true_nwk = open(f"{rd}/true_net.nwk").readline()
    troot, _ = nm.parse_enewick(true_nwk)
    true_cl = nm.softwired_clusters(troot)
    r_true = nm.num_reticulations(troot)
    gtag = {"g_500": "fasttree", "iqtree_500": "iqtree", "g_true": "true"}[gt]

    bases = {}

    def base(name):
        if name in bases:
            return bases[name]
        if name == "astral_pub":
            pf = f"{PUB}/{cond}/{rep}/astral-{gtag}.nwk"
            if os.path.exists(pf):
                nw = open(pf).readline()
            else:  # rerun ASTER's ASTRAL (identical output on the cases checked)
                subprocess.run([f"{BIN}/astral4", "-t", "4", "-i", genes, "-o", f"{work}/astral4.nwk"],
                               capture_output=True, check=True)
                nw = open(f"{work}/astral4.nwk").readline()
        elif name == "tqmc":
            subprocess.run([f"{BIN}/tree-qmc", "-i", genes, "-o", f"{work}/tqmc.nwk"], capture_output=True, check=True)
            nw = open(f"{work}/tqmc.nwk").readline()
        elif name == "wastral":
            subprocess.run([f"{BIN}/wastral", "-t", "4", "-i", genes, "-o", f"{work}/wastral.nwk"],
                           capture_output=True, check=True)
            nw = open(f"{work}/wastral.nwk").readline()
        elif name == "true_major":
            nw = nm.major_tree_newick(troot)
        elif name == "swap":
            n1 = run_variant_cached("default")["net_k1"]
            r, _ = nm.parse_enewick(n1)
            d = nm.displayed_tree_newicks(r)
            T = base("astral_pub")
            tc = nm.hardwired_clusters(nm.parse_enewick(T)[0])
            others = [x for x in d if nm.hardwired_clusters(nm.parse_enewick(x)[0]) != tc]
            nw = others[0] if others else T
        bases[name] = root_out(nw)
        return bases[name]

    cache = {}

    def run_variant_cached(v):
        if v in cache:
            return cache[v]
        bname, t, z, zmode = VARIANTS[v]
        T = base(bname)
        nets, dt, nq = run_camus(T, genes, t, z, zmode, f"{work}/{v}")
        rec = {"cond": cond, "rep": rep, "gt": gt, "variant": v, "base": bname, "t": t, "z": z, "zmode": zmode,
               "r_true": r_true, "time_s": round(dt, 3), "n_nontree_quartets": nq,
               "kmax": max(nets) if nets else 0}
        tcl = nm.softwired_clusters(nm.parse_enewick(T)[0])
        rec["base_fn"], rec["base_fp"] = nm.cluster_error(true_cl, tcl)
        for label, k in (("k1", 1), ("ktrue", r_true)):
            kk = min(k, rec["kmax"])
            nw = nets[kk][0] if kk > 0 else T
            fn, fp = nm.cluster_error(true_cl, nm.softwired_clusters(nm.parse_enewick(nw)[0]))
            rec[f"fn_{label}"], rec[f"fp_{label}"] = fn, fp
            rec[f"net_{label}"] = nw
        cache[v] = rec
        return rec

    # quartet objective (default filter) for base-tree selection
    qt = C = M = None
    with open(out, "a") as fo:
        for v in want:
            if v in done:
                continue
            rec = run_variant_cached(v)
            if VARIANTS[v][1:] == (0.5, 0, ""):
                if qt is None:
                    taxa = sorted(nm.leaves(troot))
                    qt = qq.QuartetTable(taxa)
                    C = qt.counts([l for l in open(genes) if l.strip()])
                    M = qq.filter_mask(C, 0.5)
                rec["score_k1"] = qq.network_score(qt, C, M, nm.displayed_tree_newicks(nm.parse_enewick(rec["net_k1"])[0]))
            fo.write(json.dumps(rec) + "\n")
            fo.flush()


if __name__ == "__main__":
    main()
