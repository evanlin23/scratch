"""Experiments A (impact) and B (filters) on one instance for a list of aligners.

    python run.py INSTANCE_DIR OUT.jsonl WORK_DIR --aligners mafft,magus,pasta [--detectors ts,pd,hmm]

One JSON row per (instance, aligner, condition); restartable (finished rows are skipped,
finished alignments are reused from WORK_DIR). Every score is restricted to the non-rogue
taxa C (alignment: FastSP of the induced sub-alignments; tree: FN/FP vs the true tree on C).

Conditions per aligner X:
  all          X on all taxa                     tree = FastTree on it
  all-treeC    X on all taxa, rogue rows dropped before FastTree (alignment effect only)
  oracle       X on C only (known rogues removed) tree on C
  oracle+add   oracle alignment + rogues added by HMM (UPP-style), tree on all
  D            detector D flags F; X on all-F; F added back by HMM; tree on all
  D-drop       as D but F stays out of the tree; scored on C-F (paired with all on C-F: all-on-C-F)
Aligner-independent rows: true / true-C (FastTree on the true alignment with / without rogues),
detect:D (flagged set vs known rogues).
"""

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tools  # noqa: E402
import trees  # noqa: E402
from tools import fasta  # noqa: E402

LOCK = threading.Lock()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("inst")
    p.add_argument("out")
    p.add_argument("work")
    p.add_argument("--aligners", default="mafft,magus,pasta")
    p.add_argument("--detectors", default="ts,pd,hmm")
    p.add_argument("--tree-workers", type=int, default=3)
    p.add_argument("--use-detectors", default=None, help="subset of detectors to run per aligner (default: all)")
    a = p.parse_args()
    name = os.path.basename(os.path.normpath(a.inst))
    W = os.path.join(a.work, name)
    os.makedirs(W, exist_ok=True)
    done = set()
    if os.path.exists(a.out):
        for line in open(a.out):
            r = json.loads(line)
            done.add((r["instance"], r["aligner"], r["cond"]))

    true = fasta.upper(fasta.read(os.path.join(a.inst, "true.fasta")))
    ttree = trees.read(os.path.join(a.inst, "true.tre"))
    unal = fasta.ungap(true)
    R = set(open(os.path.join(a.inst, "rogues.txt")).read().split())
    S = list(true)
    C = [n for n in S if n not in R]
    info = json.load(open(os.path.join(a.inst, "info.json")))
    base = {"instance": name, "mode": info["mode"], "k": info["k"], "length": info.get("length")}

    def emit(aligner, cond, **kw):
        row = dict(base, aligner=aligner, cond=cond, **kw)
        with LOCK:
            with open(a.out, "a") as f:
                f.write(json.dumps(row) + "\n")
            done.add((name, aligner, cond))
        print(json.dumps({k: row[k] for k in row if k not in ("flagged",)}), flush=True)

    path_locks = {}

    def tree_fn(aln, taxa_score, path):
        # trees are cached by alignment content, so identical inputs (e.g. a detector that flags
        # exactly the known rogues) share one FastTree run
        key = hashlib.md5("".join(n + aln[n] for n in sorted(aln)).encode()).hexdigest()[:16]
        os.makedirs(os.path.join(W, "trees"), exist_ok=True)
        path = os.path.join(W, "trees", key + ".tre")
        with LOCK:
            lk = path_locks.setdefault(path, threading.Lock())
        with lk:
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                tmp = "%s.%d.tmp" % (path, os.getpid())
                tools.fasttree(aln, tmp)
                os.replace(tmp, path)
        fn, fp = trees.fn_fp(ttree, trees.read(path), set(taxa_score))
        return {"FN": fn, "FP": fp}

    pool = cf.ThreadPoolExecutor(a.tree_workers)
    futs = []

    def later(aligner, cond, fn, *args, extra=None):
        if (name, aligner, cond) in done:
            return
        def job():
            try:
                r = fn(*args)
                emit(aligner, cond, **r, **(extra or {}))
            except Exception as e:  # keep going; the row is retried on restart
                print("ERROR", aligner, cond, e, flush=True)
        futs.append(pool.submit(job))

    # true-alignment trees
    later("true", "true", tree_fn, true, C, os.path.join(W, "true.tre"))
    later("true", "true-C", tree_fn, fasta.restrict(true, C), C, os.path.join(W, "trueC.tre"))

    # detectors (aligner-free inputs; initial MAFFT alignment + FastTree for TreeShrink / p-dist)
    det_path = os.path.join(W, "detect.json")
    if os.path.exists(det_path):
        det = {k: set(v) for k, v in json.load(open(det_path)).items()}
    else:
        det, timing = {}, {}
        dets = a.detectors.split(",")
        t0 = time.time()
        init, _ = tools.align("mafft", os.path.join(a.inst, "unaligned.fasta"), os.path.join(W, "init"))
        t_init = time.time() - t0
        if "pd" in dets:
            t = time.time()
            det["pd"] = tools.detect_pdist(init)[0]
            timing["pd"] = t_init + time.time() - t
        if "ts" in dets:
            t = time.time()
            tp = os.path.join(W, "init", "init.tre")
            if not os.path.exists(tp):
                tools.fasttree(init, tp)
            det["ts"] = tools.detect_treeshrink(tp, os.path.join(W, "init"))
            timing["ts"] = t_init + time.time() - t
        if "hmm" in dets:
            t = time.time()
            det["hmm"] = tools.detect_hmm(unal, os.path.join(W, "hmmdet"))[0]
            timing["hmm"] = time.time() - t
        json.dump({k: sorted(v) for k, v in det.items()}, open(det_path, "w"))
        for d, F in det.items():
            if (name, "detect", d) not in done:
                tp_ = len(F & R)
                emit("detect", d, nflag=len(F), TP=tp_, precision=tp_ / len(F) if F else None,
                     recall=tp_ / len(R), wall=round(timing[d], 1), flagged=sorted(F))

    if a.use_detectors is not None:
        det = {d: F for d, F in det.items() if d in a.use_detectors.split(",")}
    for X in a.aligners.split(","):
        conds_needed = ["all", "all-treeC", "oracle", "oracle+add"] + [d for d in det] + [d + "-drop" for d in det]
        if all((name, X, c) in done for c in conds_needed):
            continue
        cache = {}

        def get_aln(taxa, tag):
            key = hashlib.md5(",".join(sorted(taxa)).encode()).hexdigest()[:10]
            if key in cache:
                return cache[key]
            d = os.path.join(W, X, key)
            f = os.path.join(d, X + ".fasta")
            meta = os.path.join(d, "meta.json")
            if os.path.exists(meta):
                aln, wall = fasta.upper(fasta.read(f)), json.load(open(meta))["wall"]
            else:
                os.makedirs(d, exist_ok=True)
                fasta.write({n: unal[n] for n in taxa}, os.path.join(d, "in.fa"))
                aln, wall = tools.align(X, os.path.join(d, "in.fa"), d)
                json.dump({"wall": wall, "tag": tag, "n": len(taxa)}, open(meta, "w"))
            cache[key] = (aln, wall, d)
            return cache[key]

        def score_aln(aln, wall, d, tree_aln, tree_path, tree_taxa=C, extra=None):
            r = tools.fastsp_restricted(true, aln, C, os.path.join(d, "sp"))
            r.update(tree_fn(tree_aln, tree_taxa, tree_path))
            r["wall"] = wall
            r.update(extra or {})
            return r

        aln_all, w_all, d_all = get_aln(S, "all")
        later(X, "all", score_aln, aln_all, w_all, d_all, aln_all, os.path.join(d_all, "all.tre"))
        later(X, "all-treeC", score_aln, aln_all, w_all, d_all, tools.drop_allgap({n: aln_all[n] for n in C}),
              os.path.join(d_all, "allC.tre"))
        aln_o, w_o, d_o = get_aln(C, "oracle")
        later(X, "oracle", score_aln, aln_o, w_o, d_o, aln_o, os.path.join(d_o, "oracle.tre"))

        def added(aln, F, d):
            path = os.path.join(d, "added_%s.fasta" % hashlib.md5(",".join(sorted(F)).encode()).hexdigest()[:8])
            if os.path.exists(path):
                return fasta.read(path)
            t = time.time()
            m = tools.hmm_add(aln, {n: unal[n] for n in F}, os.path.join(d, "add"))
            fasta.write(m, path)
            json.dump({"wall": round(time.time() - t, 1)}, open(path + ".json", "w"))
            return m

        def add_wall(F, d):
            path = os.path.join(d, "added_%s.fasta.json" % hashlib.md5(",".join(sorted(F)).encode()).hexdigest()[:8])
            return json.load(open(path))["wall"] if os.path.exists(path) else 0.0

        if (name, X, "oracle+add") not in done:
            m = added(aln_o, sorted(R), d_o)
            later(X, "oracle+add", score_aln, m, w_o + add_wall(sorted(R), d_o), d_o, m,
                  os.path.join(d_o, "oracle_add.tre"))
        for D, F in det.items():
            F = set(F)
            ext = {"nflag": len(F), "TP": len(F & R)}
            keep = [n for n in S if n not in F]
            if not F:
                later(X, D, score_aln, aln_all, w_all, d_all, aln_all, os.path.join(d_all, "all.tre"), C, ext)
                later(X, D + "-drop", lambda: dict(tree_fn(aln_all, C, os.path.join(d_all, "all.tre")), base_FN=None, **ext))
                continue
            aln_f, w_f, d_f = get_aln(keep, D)
            if (name, X, D) not in done:
                m = added(aln_f, sorted(F), d_f)
                later(X, D, score_aln, m, w_f + add_wall(sorted(F), d_f), d_f, m, os.path.join(d_f, "added.tre"), C, ext)
            CF = [n for n in C if n not in F]

            def drop(aln_f=aln_f, d_f=d_f, CF=CF, ext=ext):
                r = tree_fn(aln_f, CF, os.path.join(d_f, "dropped.tre"))
                b = tree_fn(aln_all, CF, os.path.join(d_all, "all.tre"))
                return dict(r, base_FN=b["FN"], base_FP=b["FP"], nscored=len(CF), **ext)
            later(X, D + "-drop", drop)
    for f in futs:
        f.result()
    pool.shutdown()


if __name__ == "__main__":
    main()
