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
import anytime  # noqa: E402
import runtrees as rt  # noqa: E402
import treeerr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "runs.jsonl")
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

    def rx_cpu(self):
        """CPU seconds of this replicate's full base_raxmlng run (None if not yet run)."""
        for path in (OUT,):
            if os.path.exists(path):
                for line in open(path):
                    r = json.loads(line)
                    if (r["dataset"], r["rep"], r["aln"], r["arm"]) == (self.ds, int(self.rep), self.aln, "base_raxmlng"):
                        return r["cpu_s"]
        return None

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
                # default IQ-TREE 3, anytime (checkpoint every >= 10 s), CPU capped at this rep's RAxML-NG CPU
                cap = self.rx_cpu()
                if cap is None:
                    raise RuntimeError("base_raxmlng not done yet for this replicate (needed for the IQ-TREE cap)")
                snap = os.path.join(self.d, "trees", "anytime", self.aln + ".iqtree")
                shutil.rmtree(snap, ignore_errors=True)
                c, killed, rc = anytime.run_anytime(
                    [IQ, "-s", self.full, "-m", "GTR+G", "-T", "1", "--seed", "1", "--prefix",
                     os.path.join(work, "iq"), "-redo", "--quiet", "-cptime", "10"],
                    os.path.join(work, "iq.out"), anytime.iq_snapshot(os.path.join(work, "iq")), snap, cap_cpu=cap)
                c = Cost(c.cpu, c.wall, c.rss)
                info["killed"] = killed
                if not killed:
                    shutil.copy(os.path.join(work, "iq.treefile"), out)
                else:
                    shutil.copy(anytime.at(snap, 1e18), out)
                with open(os.path.join(snap, "index.tsv"), "a") as f:
                    f.write("%.1f\t%.1f\tfinal.tre\n" % (c.cpu, c.wall))
                shutil.copy(out, os.path.join(snap, "final.tre"))
            elif m.startswith("raxmlng"):
                # default RAxML-NG search (1 parsimony start), anytime; restartable from its checkpoint;
                # base_raxmlngcap<S>: stop once CPU > S seconds (best tree so far is the result)
                cap = float(m[len("raxmlngcap"):]) if m.startswith("raxmlngcap") else None
                snap = os.path.join(self.d, "trees", "anytime", self.aln + "." + m)
                rw = os.path.join(self.d, "trees", "rxwork_" + self.aln + "_" + m)
                os.makedirs(rw, exist_ok=True)
                pre = os.path.join(rw, "rx")
                idx = os.path.join(snap, "index.tsv")
                off = 0.0
                prev = os.path.join(rw, "spent.json")  # CPU/wall of earlier (interrupted) launches
                spent = json.load(open(prev)) if os.path.exists(prev) else {"cpu": 0.0, "wall": 0.0, "rss": 0.0}
                off = spent["cpu"]
                if os.path.exists(idx):  # launch killed from outside: its CPU is only in the snapshot index
                    off = max([off] + [float(l.split()[0]) for l in open(idx) if "final" not in l])
                    spent["cpu"] = off
                resume = os.path.exists(pre + ".raxml.ckp")
                cmd = [RX, "--search", "--msa", self.full, "--model", "GTR+G", "--threads", "1", "--seed", "1",
                       "--prefix", pre, "--tree", "pars{1}"] + ([] if resume else ["--redo"])
                spent_cost = Cost(spent["cpu"], spent["wall"], spent["rss"])
                # record progress of this launch if it gets killed from outside
                c, killed, rc = anytime.run_anytime(cmd, os.path.join(rw, "rx.out"), anytime.rx_snapshot(pre),
                                                    snap, cap_cpu=cap, cpu_offset=off)
                c = spent_cost.add(Cost(c.cpu, c.wall, c.rss))
                json.dump(c.d(), open(prev, "w"))
                info["killed"] = killed
                if not killed and rc == 0:
                    shutil.copy(pre + ".raxml.bestTree", out)
                    shutil.copy(pre + ".raxml.bestModel",
                                os.path.join(self.d, "trees", self.aln + ".base_raxmlng.model"))
                elif killed:
                    shutil.copy(anytime.at(snap, 1e18), out)
                else:
                    raise RuntimeError("raxml-ng failed, see %s" % rw)
                with open(idx, "a") as f:
                    f.write("%.1f\t%.1f\tfinal.tre\n" % (c.cpu, c.wall))
                shutil.copy(out, os.path.join(snap, "final.tre"))
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
    global OUT
    OUT = a.out
    R = Rep(a.dataset, a.rep, a.aln)
    skip = os.path.join(HERE, "skip.txt")  # "DATASET REP ARM" lines deprioritised mid-queue (compute budget)
    skip = {tuple(l.split()) for l in open(skip)} if os.path.exists(skip) else set()
    for arm in a.arms:
        if (a.dataset, int(a.rep), a.aln, arm) in done(a.out):
            continue
        if (a.dataset, str(a.rep), arm) in skip:
            print("SKIPPED (skip.txt)", a.dataset, a.rep, arm, flush=True)
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
