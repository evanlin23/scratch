"""Run species-tree methods on one gene-tree set and score them against the true species tree.

Usage: python bench.py --genes G.trees --true S.tree --key '{"data":..,"cond":..,"rep":..,"ngen":..}'
                       --out OUT.jsonl [--ngen N] [--methods m1,m2,...] [--label-mode simphy|species]
One JSON line per method: key fields + method, FN, FP, nI (internal edges of true tree), FNrate, sec.
Restartable: (key, method) pairs already in OUT are skipped. Every method runs on 1 thread.
Gene-tree leaves: 'simphy' = species before the first '_'; 'species' = leaf label is the species.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gdl", "code"))
from phylo import parse_newick, rf_error  # noqa: E402

B = "/opt/mm/root/envs/gdl/bin"
PY = B + "/python"
APRO = os.path.join(HERE, "apro")
FASTME = B + "/fastme"
ASTRALPRO = B + "/astral-pro3"
ASTRAL4 = B + "/astral4"
ASTEROID = "/opt/src/Asteroid/build/bin/asteroid"
DISCO = "/opt/src/DISCO/disco.py"
FMRFS = "/opt/src/run_fastmulrfs.sh"
WQFM = "/opt/src/run_wqfm_gdl.sh"
DUPLOSS = "/opt/src/DupLoss-2/Executables/DupLoss-2.linux"
ALL = ["astrid-multi", "astrid-pro", "astrid-pro-r0", "astrid-pro-s", "astrid-disco", "disco-astral", "astral-pro3",
       "asteroid", "fastmulrfs", "wqfm-gdl", "duploss2"]
TIMEOUT = int(os.environ.get("BENCH_TIMEOUT", "3600"))


def sh(cmd, **kw):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=TIMEOUT, **kw)


def fastme(phy, out):
    sh([FASTME, "-i", phy, "-o", out, "-m", "B", "-n", "B", "-s", "-T", "1"])


def species_relabel(genes, out, mode):
    """Gene trees with leaves renamed to bare species names (MUL-trees)."""
    import re
    with open(genes) as f, open(out, "w") as g:
        for line in f:
            if ";" not in line:
                continue
            if mode == "simphy":
                line = re.sub(r"([(,])([^(),:;]+?)_[^(),:;]*", r"\1\2", line)
            g.write(line.strip() + "\n")


def mapping(genes, out, mode):
    import re
    labs = set()
    for line in open(genes):
        for m in re.finditer(r"[(,]([^(),:;\[]+)", line):
            labs.add(m.group(1).strip())
    with open(out, "w") as f:
        for l in sorted(labs):
            f.write("%s %s\n" % (l, l.split("_")[0] if mode == "simphy" else l))


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


def run_method(m, genes, td, mode):
    u = ["-u"] if mode == "simphy" else []
    out = os.path.join(td, m + ".tre")
    if m in ("astrid-multi", "astrid-pro"):
        phy = os.path.join(td, m + ".phy")
        sh([APRO, "-i", genes, "-o", phy, "-M", m.split("-")[1]] + u)
        fastme(phy, out)
    elif m == "astrid-pro-r0":  # ablation: gene-tree root not counted
        phy = os.path.join(td, m + ".phy")
        sh([APRO, "-i", genes, "-o", phy, "-M", "pro", "-R", "0"] + u)
        fastme(phy, out)
    elif m == "astrid-pro-s":
        phy1, t1, phy = (os.path.join(td, x) for x in ("s1.phy", "s1.tre", "s.phy"))
        sh([APRO, "-i", genes, "-o", phy1, "-M", "pro"] + u)
        fastme(phy1, t1)
        sh([APRO, "-i", genes, "-o", phy, "-M", "pros", "-s", t1] + u)
        fastme(phy, out)
    elif m in ("astrid-disco", "disco-astral"):
        dec = os.path.join(td, "disco.trees")
        if not os.path.exists(dec):
            cmd = [PY, DISCO, "-i", genes, "-o", dec]
            if mode == "simphy":
                cmd += ["-d", "_"]
            sh(cmd)
            strip_root_leaf(dec)
        if m == "astrid-disco":
            phy = os.path.join(td, "disco.phy")
            sh([APRO, "-i", dec, "-o", phy, "-M", "multi"])  # single-copy trees: ASTRID
            fastme(phy, out)
        else:
            sh([ASTRAL4, "-i", dec, "-o", out, "-t", "1"])
    elif m == "astral-pro3":
        mp = os.path.join(td, "map.txt")
        mapping(genes, mp, mode)
        sh([ASTRALPRO, "-i", genes, "-a", mp, "-o", out, "-t", "1", "-u", "0"])
    elif m == "asteroid":
        mp = os.path.join(td, "map.txt")
        mapping(genes, mp, mode)
        sh([ASTEROID, "-i", genes, "-m", mp, "-p", os.path.join(td, "ast")])
        shutil.copy(os.path.join(td, "ast.bestTree.newick"), out)
    elif m == "fastmulrfs":
        if mode == "simphy":
            sh(["bash", FMRFS, genes, out])
        else:
            mp = os.path.join(td, "map.txt")
            mapping(genes, mp, mode)
            sh(["bash", FMRFS, genes, mp, out])
    elif m == "wqfm-gdl":
        mp = "-"
        if mode != "simphy":
            mp = os.path.join(td, "map.txt")
            mapping(genes, mp, mode)
        sh(["bash", WQFM, genes, out, mp], env=dict(os.environ, WQFM_MEM="6g"))
    elif m == "duploss2":
        sp = os.path.join(td, "sp.trees")
        species_relabel(genes, sp, mode)
        sh([DUPLOSS, "-i", sp, "-o", out])
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
    ap.add_argument("--ngen", type=int, default=0)
    ap.add_argument("--methods", default=",".join(ALL))
    ap.add_argument("--label-mode", default="simphy")
    a = ap.parse_args()
    key = json.loads(a.key)
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            r = json.loads(l)
            done.add((json.dumps({k: r[k] for k in key}, sort_keys=True), r["method"]))
    ks = json.dumps(key, sort_keys=True)
    todo = [m for m in a.methods.split(",") if (ks, m) not in done]
    if not todo:
        return
    true = parse_newick(open(a.true).read().strip().split("\n")[0])
    with tempfile.TemporaryDirectory(dir=os.environ.get("BENCH_TMP")) as td:
        genes = os.path.join(td, "genes.trees")
        with open(a.genes) as f, open(genes, "w") as g:
            n = 0
            for line in f:
                if ";" not in line:
                    continue
                g.write(line.strip() + "\n")
                n += 1
                if a.ngen and n >= a.ngen:
                    break
        for m in todo:
            t0 = time.time()
            rec = dict(key, method=m)
            try:
                nwk = run_method(m, genes, td, a.label_mode)
                rec["sec"] = round(time.time() - t0, 3)
                est = parse_newick(nwk)
                fn, fp, i1, i2 = rf_error(est, true)
                rec.update(FN=fn, FP=fp, nI=i1, FNrate=round(fn / i1, 5))
            except Exception as e:  # record failures (timeouts, crashes) so reruns skip them
                rec["sec"] = round(time.time() - t0, 3)
                rec["error"] = type(e).__name__
            with open(a.out, "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
