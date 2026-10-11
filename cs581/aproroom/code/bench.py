"""Run species-tree methods on one gene-tree set and score them against the true species tree.
(aproroom version: adds true-tag variants (*-tt, need --tagged), thread control (--threads, for astral-pro3 and
wqfm-gdl's JVM), the ASTRAL-Pro3 hybrids (apro3-guide*, apro3-fast) and peak RSS per method (os.wait4).)

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
sys.path.insert(0, HERE)
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
ALL = ["astrid-multi", "astrid-pro", "astrid-disco", "disco-astral", "astral-pro3", "asteroid", "wqfm-gdl",
       "astrid-pro-tt", "astrid-disco-tt", "disco-astral-tt"]
DISCO_TT = os.path.join(HERE, "disco_tt.py")
THREADS = 1
PEAK = [0]
TIMEOUT = int(os.environ.get("BENCH_TIMEOUT", "3600"))


def sh(cmd, **kw):
    """run cmd; record its peak RSS (KB, incl. waited-for descendants) in PEAK[0]; kill the group on timeout"""
    import signal
    p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True, **kw)
    t0 = time.time()
    while True:
        pid, st, ru = os.wait4(p.pid, os.WNOHANG)
        if pid:
            break
        if time.time() - t0 > TIMEOUT:
            os.killpg(p.pid, signal.SIGKILL)
            os.wait4(p.pid, 0)
            raise subprocess.TimeoutExpired(cmd, TIMEOUT)
        time.sleep(0.05)
    p.returncode = os.waitstatus_to_exitcode(st)
    PEAK[0] = max(PEAK[0], ru.ru_maxrss)
    if p.returncode:
        raise subprocess.CalledProcessError(p.returncode, cmd)


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


def rf_shared(est, true):
    """(FN, FP, nI_true, nI_est) on the leaf set shared by both trees (empirical references may be
    induced on fewer species); identical to phylo.rf_error when the leaf sets agree."""
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
    return len(bt - be), len(be - bt), len(bt), len(be)


def run_method(m, genes, td, mode, tagged=None):
    u = ["-u"] if mode == "simphy" else []
    out = os.path.join(td, m + ".tre")
    if m in ("astrid-multi", "astrid-pro"):
        phy = os.path.join(td, m + ".phy")
        sh([APRO, "-i", genes, "-o", phy, "-M", m.split("-")[1], "-t", str(THREADS)] + u)
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
    elif m == "astrid-pro-tt":  # true root and tags
        phy = os.path.join(td, m + ".phy")
        sh([APRO, "-i", tagged, "-o", phy, "-M", "pro", "-T"] + u)
        fastme(phy, out)
    elif m in ("astrid-disco-tt", "disco-astral-tt"):  # DISCO decomposition at the TRUE duplications
        dec = os.path.join(td, "disco_tt.trees")
        if not os.path.exists(dec):
            sh([PY, DISCO_TT, tagged, dec])
        if m == "astrid-disco-tt":
            phy = os.path.join(td, "disco_tt.phy")
            sh([APRO, "-i", dec, "-o", phy, "-M", "multi"])
            fastme(phy, out)
        else:
            sh([ASTRAL4, "-i", dec, "-o", out, "-t", str(THREADS)])
    elif m in ("astral-pro3", "apro3-fast", "apro3-guide", "apro3-guide-fast"):
        mp = os.path.join(td, "map.txt")
        mapping(genes, mp, mode)
        extra = []
        if m.startswith("apro3-guide"):
            g = os.path.join(td, "apro_guide.tre")
            if not os.path.exists(g):
                phy = os.path.join(td, "apro_guide.phy")
                sh([APRO, "-i", genes, "-o", phy, "-M", "pro"] + u)
                fastme(phy, g)
            extra += ["-g", g]
        if m.endswith("fast"):
            extra += ["-r", "1", "-s", "0"]
        sh([ASTRALPRO, "-i", genes, "-a", mp, "-o", out, "-t", str(THREADS), "-u", "0"] + extra)
    elif m == "asteroid":
        mp = []
        if mode == "simphy":  # with species labels Asteroid needs no mapping (it rejects identity maps)
            mp = ["-m", os.path.join(td, "map.txt")]
            mapping(genes, mp[1], mode)
        g4 = os.path.join(td, "genes4.trees")  # Asteroid cannot parse gene trees with < 4 leaves
        with open(genes) as f, open(g4, "w") as g:
            for line in f:
                if line.count(",") >= 3:
                    g.write(line)
        sh([ASTEROID, "-i", g4, "-p", os.path.join(td, "ast")] + mp)
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
        sh(["bash", WQFM, genes, out], env=dict(os.environ, WQFM_MEM=os.environ.get("WQFM_MEM", "4g"),
                                               JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=%d" % THREADS))
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
    ap.add_argument("--tagged", default=None, help="true-tagged gene trees (truetag.py), same genes/order")
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    global THREADS
    THREADS = a.threads
    key = json.loads(a.key)
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            r = json.loads(l)
            done.add((json.dumps({k: r.get(k) for k in key}, sort_keys=True), r["method"]))
    ks = json.dumps(key, sort_keys=True)
    todo = [m for m in a.methods.split(",") if (ks, m) not in done]
    if not todo:
        return
    txt = "".join(l.strip() for l in open(a.true))
    true = parse_newick(txt[:txt.index(";") + 1].replace("[&R]", "").replace("'", ""))
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
        tagged = None
        if a.tagged:
            tagged = os.path.join(td, "tagged.trees")
            with open(a.tagged) as f, open(tagged, "w") as g:
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
            rec = dict(key, method=m, threads=THREADS)
            PEAK[0] = 0
            try:
                if m.endswith("-tt") and not tagged:
                    continue
                nwk = run_method(m, genes, td, a.label_mode, tagged)
                rec["rssMB"] = round(PEAK[0] / 1024, 1)
                rec["sec"] = round(time.time() - t0, 3)
                est = parse_newick(nwk)
                fn, fp, i1, i2 = rf_shared(est, true)
                rec.update(FN=fn, FP=fp, nI=i1, FNrate=round(fn / i1, 5))
            except Exception as e:  # record failures (timeouts, crashes) so reruns skip them
                rec["sec"] = round(time.time() - t0, 3)
                rec["rssMB"] = round(PEAK[0] / 1024, 1)
                rec["error"] = type(e).__name__
            with open(a.out, "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
