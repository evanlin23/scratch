"""Backbone -> place -> graft -> polish pipelines, constrained ML, and baselines, on one replicate.

    python pipe.py DATASET REP ARM [ARM ...] [--aln true_align] [--out ../results/runs.jsonl]

ARM names:
  base_fasttree | base_iqfast | base_iqtree | base_raxmlng      (all sequences)
  constr_<BB>_<TAU>                       (A) RAxML-NG (GTR+G, 1 parsimony start) on all sequences with
                                          the backbone tree as a non-comprehensive --tree-constraint
  place_<BB>_<TAU>_<EPA>_<POLISH>         (B/C) backbone tree -> RAxML-NG --evaluate (branch lengths +
                                          GTR+G for placement) -> EPA-ng (EPA = fix: patched, stock:
                                          v0.3.8) places every non-backbone sequence -> gappa examine
                                          graft (best-LWR placement, --fully-resolve) -> POLISH:
                                            graft  : no polish (the grafted tree itself)
                                            rxfast : RAxML-NG from the grafted tree, fast mode
                                                     (--opt-topology simplified --stop-rule kh-mult)
                                            rxfull : RAxML-NG from the grafted tree, default search
                                            iqfast : IQ-TREE 3 -t grafted --fast
                                            ft     : FastTree -intree grafted (NNI + SPR rounds)
  BB = ft (FastTree -gtr -gamma) | iqf (IQ-TREE 3 --fast); TAU in {0.5, 0.75}: backbone = sequences with
  ungapped length >= TAU * median ungapped length (median over all sequences).

Every tool runs with 1 thread. Shared steps (backbone tree, placement + graft) are cached per replicate
under trees/cache/ with their measured cost; every arm's cpu/wall include the cost of all of its steps,
as if run from scratch. peak_rss_mb = max RSS over the arm's steps. A row is appended per finished arm;
arms already in --out are skipped (restartable).
"""
import argparse
import fcntl
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402
import treeerr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FT = "/usr/bin/FastTree"
IQ = "/opt/mm/root/envs/bio/bin/iqtree3"
RX = "/opt/mm/root/envs/fml/bin/raxml-ng"
GAPPA = "/opt/mm/root/envs/fml/bin/gappa"
EPA = {"fix": "/opt/src/epa-ng-fix/bin/epa-ng", "stock": "/opt/src/epa-ng/bin/epa-ng"}


class Cost:
    def __init__(self, cpu=0.0, wall=0.0, rss=0.0):
        self.cpu, self.wall, self.rss = cpu, wall, rss

    def add(self, o):
        return Cost(self.cpu + o.cpu, self.wall + o.wall, max(self.rss, o.rss))

    def d(self):
        return {"cpu": round(self.cpu, 1), "wall": round(self.wall, 1), "rss": round(self.rss, 1)}


def run(cmd, log, stdout=None, cwd=None):
    t = time.time()
    with open(log, "w") as lf:
        p = subprocess.Popen(cmd, stdout=stdout or lf, stderr=lf, cwd=cwd)
        _, status, ru = os.wait4(p.pid, 0)
    if os.waitstatus_to_exitcode(status) != 0:
        raise RuntimeError("failed: %s (log %s)" % (" ".join(cmd), log))
    return Cost(ru.ru_utime + ru.ru_stime, time.time() - t, ru.ru_maxrss / 1024.0)


def fasttree(aln, out, work, extra=()):
    with open(out, "w") as o:
        return run([FT, "-nt", "-gtr", "-gamma"] + list(extra) + [aln], os.path.join(work, "ft.log"), stdout=o)


def iqtree(aln, out, work, extra=()):
    c = run([IQ, "-s", aln, "-m", "GTR+G", "-T", "1", "--seed", "1", "--prefix", os.path.join(work, "iq"),
             "-redo", "--quiet"] + list(extra), os.path.join(work, "iq.out"))
    shutil.copy(os.path.join(work, "iq.treefile"), out)
    return c


def raxml(aln, out, work, extra=(), model_out=None):
    c = run([RX, "--search", "--msa", aln, "--model", "GTR+G", "--threads", "1", "--seed", "1", "--prefix",
             os.path.join(work, "rx"), "--redo"] + list(extra), os.path.join(work, "rx.out"))
    shutil.copy(os.path.join(work, "rx.raxml.bestTree"), out)
    if model_out:
        shutil.copy(os.path.join(work, "rx.raxml.bestModel"), model_out)
    return c


class Rep:
    def __init__(self, ds, rep, aln):
        self.ds, self.rep, self.aln = ds, rep, aln
        self.d = os.path.join(rt.MLDATA, ds, "R%s" % rep)
        self.true = os.path.join(self.d, "true_tree.tre")
        self.cache = os.path.join(self.d, "trees", "cache")
        os.makedirs(self.cache, exist_ok=True)
        self.full = os.path.join(self.d, aln + ".clean.fasta")
        with self.lock("clean"):
            if not os.path.exists(self.full):
                rt.clean_alignment(os.path.join(self.d, aln + ".fasta"), self.full + ".tmp")
                os.replace(self.full + ".tmp", self.full)
        self.names, self.seqs = rt.read_fasta(self.full)
        self.lens = [len(s.replace("-", "")) for s in self.seqs]
        self.med = statistics.median(self.lens)

    def lock(self, key):
        class L:
            def __enter__(s):
                s.f = open(os.path.join(self.cache, key + ".lock"), "w")
                fcntl.flock(s.f, fcntl.LOCK_EX)

            def __exit__(s, *a):
                fcntl.flock(s.f, fcntl.LOCK_UN)
                s.f.close()
        return L()

    def cached(self, key, fn):
        """fn(work, outprefix) -> Cost; result files <cache>/<key>.*; cost stored in <key>.cost.json"""
        with self.lock(key):
            cj = os.path.join(self.cache, key + ".cost.json")
            if os.path.exists(cj):
                c = json.load(open(cj))
                return Cost(c["cpu"], c["wall"], c["rss"])
            work = tempfile.mkdtemp(prefix="fml_%s_%s_%s_" % (self.ds, self.rep, key))
            c = fn(work, os.path.join(self.cache, key))
            json.dump(c.d(), open(cj + ".tmp", "w"))
            os.replace(cj + ".tmp", cj)
            shutil.rmtree(work, ignore_errors=True)
            return c

    def split(self, tau):
        bb = [i for i, l in enumerate(self.lens) if l >= tau * self.med]
        return bb, [i for i in range(len(self.names)) if i not in set(bb)]

    def backbone(self, bbm, tau):
        key = "%s.bb_%s_%g" % (self.aln, bbm, tau)

        def fn(work, pre):
            bb, _ = self.split(tau)
            rt.write_fasta(pre + ".raw.fasta", [self.names[i] for i in bb], [self.seqs[i] for i in bb])
            rt.clean_alignment(pre + ".raw.fasta", pre + ".fasta")
            if bbm == "ft":
                return fasttree(pre + ".fasta", pre + ".tre", work)
            if bbm == "iqf":
                return iqtree(pre + ".fasta", pre + ".tre", work, ["--fast"])
            raise ValueError(bbm)
        return key, self.cached(key, fn)

    def placed(self, bbm, tau, epa):
        bkey, bcost = self.backbone(bbm, tau)
        bpre = os.path.join(self.cache, bkey)
        # branch lengths + GTR+G parameters of the backbone tree for EPA-ng (shared by fix/stock)
        ekey = bkey + ".eval"

        def fe(work, pre):
            c = run([RX, "--evaluate", "--msa", bpre + ".fasta", "--tree", bpre + ".tre", "--model", "GTR+G",
                     "--threads", "1", "--seed", "1", "--prefix", os.path.join(work, "ev"), "--redo"],
                    os.path.join(work, "ev.out"))
            shutil.copy(os.path.join(work, "ev.raxml.bestTree"), pre + ".tre")
            shutil.copy(os.path.join(work, "ev.raxml.bestModel"), pre + ".model")
            return c
        ecost = self.cached(ekey, fe)
        epre = os.path.join(self.cache, ekey)
        pkey = "%s.place_%s" % (bkey, epa)

        def fp(work, pre):
            bb, fr = self.split(tau)
            ref, qry = os.path.join(work, "ref.fasta"), os.path.join(work, "qry.fasta")
            rt.write_fasta(ref, [self.names[i] for i in bb], [self.seqs[i] for i in bb])
            rt.write_fasta(qry, [self.names[i] for i in fr], [self.seqs[i] for i in fr])
            c = run([EPA[epa], "--ref-msa", ref, "--tree", epre + ".tre", "--query", qry, "--model",
                     epre + ".model", "--threads", "1", "--outdir", work, "--redo"], os.path.join(work, "epa.log"))
            shutil.copy(os.path.join(work, "epa_result.jplace"), pre + ".jplace")
            c = c.add(run([GAPPA, "examine", "graft", "--jplace-path", pre + ".jplace", "--fully-resolve",
                           "--out-dir", work, "--allow-file-overwriting"], os.path.join(work, "graft.log")))
            shutil.copy(os.path.join(work, os.path.basename(pre) + ".newick"), pre + ".tre")
            return c
        pcost = self.cached(pkey, fp)
        return os.path.join(self.cache, pkey) + ".tre", bcost.add(ecost).add(pcost), bpre + ".tre"

    def arm(self, name):
        out = os.path.join(self.d, "trees", "%s.%s.tre" % (self.aln, name))
        work = tempfile.mkdtemp(prefix="fml_%s_%s_%s_" % (self.ds, self.rep, name))
        p = name.split("_")
        info = {}
        if p[0] == "base":
            m = "_".join(p[1:])
            if m == "fasttree":
                c = fasttree(self.full, out, work)
            elif m == "iqfast":
                c = iqtree(self.full, out, work, ["--fast"])
            elif m == "iqtree":
                c = iqtree(self.full, out, work)
            elif m == "raxmlng":
                c = raxml(self.full, out, work, ["--tree", "pars{1}"],
                          model_out=os.path.join(self.d, "trees", self.aln + ".base_raxmlng.model"))
            else:
                raise ValueError(name)
        elif p[0] == "constr":
            bbm, tau = p[1], float(p[2])
            bkey, c = self.backbone(bbm, tau)
            tb = os.path.join(self.cache, bkey + ".tre")
            info["backbone_fn"] = treeerr.error(self.true, tb)["fn_rate"]
            c = c.add(raxml(self.full, out, work, ["--tree", "pars{1}", "--tree-constraint", tb]))
        elif p[0] == "place":
            bbm, tau, epa, pol = p[1], float(p[2]), p[3], p[4]
            g, c, tb = self.placed(bbm, tau, epa)
            info["backbone_fn"] = treeerr.error(self.true, tb)["fn_rate"]
            info["graft_fn"] = treeerr.error(self.true, g)["fn_rate"]
            if pol == "graft":
                shutil.copy(g, out)
            elif pol == "rxfast":
                c = c.add(raxml(self.full, out, work, ["--tree", g, "--opt-topology", "simplified",
                                                       "--stop-rule", "kh-mult"]))
            elif pol == "rxfull":
                c = c.add(raxml(self.full, out, work, ["--tree", g]))
            elif pol == "iqfast":
                c = c.add(iqtree(self.full, out, work, ["-t", g, "--fast"]))
            elif pol == "ft":
                c = c.add(fasttree(self.full, out, work, ["-intree", g]))
            else:
                raise ValueError(name)
        else:
            raise ValueError(name)
        shutil.rmtree(work, ignore_errors=True)
        e = treeerr.error(self.true, out)
        return {"dataset": self.ds, "rep": int(self.rep), "aln": self.aln, "arm": name,
                "fn": e["fn_rate"], "fp": e["fp_rate"], "rf": e["rf_rate"], "est_int": e["est_int"],
                "cpu_s": round(c.cpu, 1), "wall_s": round(c.wall, 1), "peak_rss_mb": round(c.rss, 1),
                "tree": out, **info}


def done(path):
    s = set()
    if os.path.exists(path):
        for line in open(path):
            r = json.loads(line)
            s.add((r["dataset"], r["rep"], r["aln"], r["arm"]))
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("rep")
    ap.add_argument("arms", nargs="+")
    ap.add_argument("--aln", default="true_align")
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "runs.jsonl"))
    a = ap.parse_args()
    R = Rep(a.dataset, a.rep, a.aln)
    for arm in a.arms:
        if (a.dataset, int(a.rep), a.aln, arm) in done(a.out):
            continue
        try:
            row = R.arm(arm)
        except Exception as ex:  # keep going with the other arms
            print("FAILED", a.dataset, a.rep, arm, ex, flush=True)
            continue
        with open(a.out, "a") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
