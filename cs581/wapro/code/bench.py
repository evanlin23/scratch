"""Run species-tree methods on one gene-tree set and score them against the true species tree.

Usage: python bench.py --genes G.trees --true S.tree --key '{"data":..,"cond":..,"rep":..}' --out OUT.jsonl
                       [--methods m1,m2,...] [--threads 1]
One JSON line per method: key fields + method, threads, FN, FP, nI, FNrate, wall (s), cpu (s, children).
Restartable: (key, method, threads) triples already in OUT are skipped.
Gene-tree leaves are SimPhy labels 'species_locus_individual' (species = text before the first '_').

Methods
  wapro:<flags>     weighted ASTRID-Pro, e.g. 'wapro:-W supe' (flags are passed to code/wapro), + FastME 2
  astrid-pro        = wapro with unit weights (the published-round ASTRID-Pro)
  astral-pro3       ASTER v1.25.3.8 (no weighting option exists for multi-copy input)
  disco-wastral     DISCO v1.4.1 then wASTRAL (ASTER v1.25.3.8, --mode 1 = hybrid weighting, support in [0,1])
  disco-wastrid     DISCO then wASTRID-s (internode v0.0.7, -m support -b 0-1)
  astrid-disco      DISCO then ASTRID (internode -m internode, unweighted)
  asteroid          Asteroid (missing-data correction on), internode distances
  asteroid-bl       Asteroid --use-gene-bl (its only length/support-aware setting)
  wqfm-gdl          wQFM-GDL v1.0.3, tree mode (-t)
"""
import argparse
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from phylo import parse_newick  # noqa: E402

B = "/opt/mm/root/envs/gdl/bin"
PY = B + "/python"
WAPRO = os.path.join(HERE, "wapro")
FASTME = B + "/fastme"
ASTRALPRO = B + "/astral-pro3"
WASTRAL = B + "/wastral"
WASTRID = "/opt/src/internode/target/release/wastrid"
ASTEROID = "/opt/src/Asteroid/build/bin/asteroid"
DISCO = "/opt/src/DISCO/disco.py"
WQFM = "/opt/src/wQFM-GDL/wQFM-GDL.sh"
TIMEOUT = int(os.environ.get("BENCH_TIMEOUT", "3600"))


def sh(cmd, threads=1, **kw):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), PATH=B + ":" + os.environ["PATH"])
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=TIMEOUT, env=env, **kw)


def fastme(phy, out, threads=1):
    sh([FASTME, "-i", phy, "-o", out, "-m", "B", "-n", "B", "-s", "-T", str(threads)], threads)


def mapping(genes, out):
    labs = set()
    for line in open(genes):
        for m in re.finditer(r"[(,]([^(),:;\[]+)", line):
            labs.add(m.group(1).strip())
    with open(out, "w") as f:
        for l in sorted(labs):
            f.write("%s %s\n" % (l, l.split("_")[0]))


def strip_root_leaf(path):
    """DISCO (v1.4.1) sometimes emits a leaf labelled ROOT; prune it."""
    import treeswift
    out = []
    for line in open(path):
        if ";" not in line:
            continue
        if "ROOT" in line:
            t = treeswift.read_tree_newick(line)
            labs = {x.label for x in t.traverse_leaves()} - {"ROOT"}
            if len(labs) < 4:
                continue
            t = t.extract_tree_with(labs)
            t.suppress_unifurcations()
            line = t.newick()
        out.append(line.strip())
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def rf(est, true):
    """(FN, FP, nI_true) on the shared leaf set."""
    shared = sorted(set(est.label[v] for v in est.leaves()) & set(true.label[v] for v in true.leaves()))
    idx = {x: i for i, x in enumerate(shared)}
    n = len(shared)
    full = (1 << n) - 1

    def bip(t):
        mask = [0] * len(t.parent)
        for v in t.postorder():
            if t.children[v]:
                for c in t.children[v]:
                    mask[v] |= mask[c]
            elif t.label[v] in idx:
                mask[v] = 1 << idx[t.label[v]]
        out = set()
        for m in mask:
            if m & 1:
                m = full ^ m
            if 2 <= bin(m).count("1") <= n - 2:
                out.add(m)
        return out
    be, bt = bip(est), bip(true)
    return len(bt - be), len(be - bt), len(bt)


def disco(genes, td):
    """Run DISCO afresh for every method that uses it, so its time is included in each."""
    dec = os.path.join(td, "disco.trees")
    sh([PY, DISCO, "-i", genes, "-o", dec, "-d", "_"])
    strip_root_leaf(dec)
    return dec


def run_method(m, genes, td, th):
    out = os.path.join(td, re.sub(r"[^A-Za-z0-9]", "_", m) + ".tre")
    if m.startswith("wapro:") or m == "astrid-pro":
        flags = m[6:].split() if m.startswith("wapro:") else []
        phy = out + ".phy"
        sh([WAPRO, "-i", genes, "-o", phy, "-M", "pro", "-u"] + flags)
        fastme(phy, out, th)
    elif m in ("disco-wastrid", "astrid-disco"):
        dec = disco(genes, td)
        mode = ["-m", "support", "-b", "0-1"] if m == "disco-wastrid" else ["-m", "internode"]
        sh([WASTRID, "-i", dec, "-o", out, "-t", str(th)] + mode, th)
    elif m == "disco-wastral":
        dec = disco(genes, td)
        sh([WASTRAL, "-i", dec, "-o", out, "-t", str(th), "--mode", "1", "-u", "0"], th)
    elif m == "astral-pro3":
        mp = os.path.join(td, "map.txt")
        mapping(genes, mp)
        sh([ASTRALPRO, "-i", genes, "-a", mp, "-o", out, "-t", str(th), "-u", "0"], th)
    elif m in ("asteroid", "asteroid-bl"):
        mp = os.path.join(td, "map.txt")
        mapping(genes, mp)
        g4 = os.path.join(td, "genes4.trees")  # Asteroid cannot parse gene trees with < 4 leaves
        with open(genes) as f, open(g4, "w") as g:
            for line in f:
                if line.count(",") >= 3:
                    g.write(line)
        extra = ["--use-gene-bl"] if m == "asteroid-bl" else []
        sh([ASTEROID, "-i", g4, "-m", mp, "-p", os.path.join(td, "ast" + m)] + extra, th)
        shutil.copy(os.path.join(td, "ast" + m + ".bestTree.newick"), out)
    elif m == "wqfm-gdl":
        # wQFM-GDL accepts 'species_copy' labels with a single '_' only
        g1 = os.path.join(td, "wq.trees")
        with open(genes) as f, open(g1, "w") as g:
            for line in f:
                g.write(re.sub(r"([(,])([^(),:;_]+)_([^(),:;_]+)_([^(),:;]+)", r"\1\2_\3c\4", line))
        # its pipeline calls ./scripts/paup, so it must run from its own directory
        sh(["bash", WQFM, "-i", g1, "-o", out, "-t", "-m", "4g"], th, cwd=os.path.dirname(WQFM))
    else:
        raise ValueError(m)
    with open(out) as f:
        return next(l.strip() for l in f if l.strip().startswith("(") and ";" in l)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True)
    ap.add_argument("--true", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--methods", required=True)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    key = json.loads(a.key)
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            r = json.loads(l)
            done.add((json.dumps({k: r[k] for k in key}, sort_keys=True), r["method"], r.get("threads", 1)))
    ks = json.dumps(key, sort_keys=True)
    true = parse_newick(open(a.true).read().strip())
    os.makedirs("/opt/tmp/bench", exist_ok=True)
    with tempfile.TemporaryDirectory(dir="/opt/tmp/bench") as td:
        for m in a.methods.split(","):
            if (ks, m, a.threads) in done:
                continue
            rec = dict(key, method=m, threads=a.threads)
            r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
            t0 = time.time()
            try:
                nwk = run_method(m, a.genes, td, a.threads)
                wall = time.time() - t0
                r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
                fn, fp, ni = rf(parse_newick(nwk), true)
                rec.update(FN=fn, FP=fp, nI=ni, FNrate=fn / ni, wall=round(wall, 3),
                           cpu=round(r1.ru_utime + r1.ru_stime - r0.ru_utime - r0.ru_stime, 3))
            except Exception as e:  # noqa: BLE001
                rec.update(error=type(e).__name__, wall=round(time.time() - t0, 1))
            with open(a.out, "a") as f:
                f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
