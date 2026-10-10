"""Isolate merger error: merge two halves of a replicate with each merger and score.

    python run.py RESULTS.jsonl DATASET:REP_DIR [DATASET:REP_DIR ...] [--jobs 4] [--conditions oracle,fftnsi]
                  [--mergers mafft-merge,muscle3,opal,gcm,progdp]

Per replicate (work dir /opt/runs/pairmerge/<dataset>/<rep>):
  * split the taxa into two halves at the centroid edge of the TRUE tree;
  * sub-alignments: `oracle` = true alignment restricted to each half (merger error only),
    `fftnsi` = MAFFT FFT-NS-i (--retree 2 --maxiterate 2) on each half
    (L-INS-i on 500 x ~1000-bp sequences takes >30 core-minutes per half: over budget);
  * every merger merges the same two sub-alignments; GCM runs first and its
    MAGUS-built backbones are reused by progdp (same evidence, only the trace differs).
Scores: FastSP SPFN/SPFP on the full alignment, plus SPFN/SPFP restricted to
cross-half residue pairs (the only pairs a merger decides). Rows already in
RESULTS.jsonl are skipped, so the script is restartable.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import mergers
import tree
from gcmx import fasta, score

RUNS = os.environ.get("PAIRMERGE_RUNS", "/opt/runs/pairmerge")
ALL_MERGERS = ["mafft-merge", "muscle3", "opal", "gcm", "progdp"]


def find_inputs(rep_dir):
    """(true alignment, true tree) of a MAGUS-paper replicate directory."""
    for aln, tre in (("rose.aln.true.fasta", "rose.tt"), ("true_align.txt", "true_tree.tre")):
        if os.path.exists(os.path.join(rep_dir, aln)):
            return os.path.join(rep_dir, aln), os.path.join(rep_dir, tre)
    raise FileNotFoundError("no true alignment/tree in " + rep_dir)


def cross_scores(true, est, A):
    """FN/FP rates over homology pairs (a, b) with a in A, b not in A."""
    taxa = list(true)
    tcol, ecol, side = [], [], []
    for t in taxa:
        ts, es = true[t], est[t]
        tc = np.array([i for i, c in enumerate(ts) if c not in "-."])
        ec = np.array([i for i, c in enumerate(es) if c not in "-."])
        assert len(tc) == len(ec), t
        tcol.append(tc); ecol.append(ec); side.append(np.full(len(tc), t in A))
    tcol, ecol, side = np.concatenate(tcol), np.concatenate(ecol), np.concatenate(side)

    def cross_pairs(*keys):
        key = np.zeros(len(side), dtype=np.int64)
        for k in keys:
            key = key * (k.max() + 1) + k
        _, inv = np.unique(key, return_inverse=True)
        na = np.bincount(inv, weights=side.astype(float))
        nb = np.bincount(inv, weights=(~side).astype(float))
        return float((na * nb).sum())

    t, e, shared = cross_pairs(tcol), cross_pairs(ecol), cross_pairs(tcol, ecol)
    return {"crossFN": 1 - shared / t, "crossFP": 1 - shared / e if e else 0.0,
            "crossTrue": t, "crossEst": e}


def prepare(dataset, rep_dir, conditions):
    work = os.path.join(RUNS, dataset)
    os.makedirs(work, exist_ok=True)
    aln_path, tree_path = find_inputs(rep_dir)
    true = fasta.upper(fasta.read(aln_path))
    fasta.write(true, os.path.join(work, "true.fasta"))
    A, B = tree.centroid_split(open(tree_path).read())
    assert set(A) | set(B) == set(true), "tree/alignment taxa differ"
    info = {"nA": len(A), "nB": len(B)}
    for name, half in (("A", A), ("B", B)):
        fasta.write(fasta.restrict(true, half), os.path.join(work, "oracle_{}.fa".format(name)))
        fasta.write(fasta.ungap({t: true[t] for t in half}), os.path.join(work, "unaligned_{}.fa".format(name)))
        out = os.path.join(work, "fftnsi_{}.fa".format(name))
        if "fftnsi" in conditions and not os.path.exists(out):
            start = time.time()
            with open(out + ".tmp", "w") as o:
                subprocess.run(["mafft", "--retree", "2", "--maxiterate", "2", "--thread", "1", "--quiet",
                                os.path.join(work, "unaligned_{}.fa".format(name))], stdout=o, check=True)
            # MAFFT writes lower case; MAGUS matches letters case-sensitively against its upper-case backbones
            fasta.write(fasta.upper(fasta.read(out + ".tmp")), out)
            os.remove(out + ".tmp")
            info["fftnsi_seconds_" + name] = round(time.time() - start, 1)
    return work, true, set(A), info


def run_replicate(results, dataset, rep_dir, conditions, merger_names, done):
    records = []
    work, true, A, info = prepare(dataset, rep_dir, conditions)
    # MAGUS backbones depend only on the unaligned sequences: build them once (first gcm run)
    # and reuse them for every condition, so gcm and progdp always see identical evidence.
    bb = os.path.join(work, "backbones")
    for cond in conditions:
        a = os.path.join(work, "{}_A.fa".format(cond))
        b = os.path.join(work, "{}_B.fa".format(cond))
        for m in merger_names:
            if (dataset, cond, m) in done:
                continue
            mwork = os.path.join(work, cond + "_" + m)
            shutil.rmtree(mwork, ignore_errors=True)
            os.makedirs(mwork)
            out = os.path.join(work, "{}_{}.fasta".format(cond, m))
            rec = {"dataset": dataset, "condition": cond, "merger": m, **info}
            start = time.time()
            try:
                if m == "mafft-merge":
                    mergers.mafft_merge(a, b, out, mwork, 1)
                elif m == "muscle3":
                    mergers.muscle3(a, b, out, mwork, 1)
                elif m == "opal":
                    mergers.opal(a, b, out, mwork, 1)
                elif m == "gcm":
                    if glob.glob(os.path.join(bb, "*")):
                        mergers.gcm(a, b, out, mwork, 1, backbones=bb)
                    else:
                        shutil.rmtree(bb + ".tmp", ignore_errors=True)
                        mergers.gcm(a, b, out, mwork, 1, keep_backbones=bb + ".tmp")
                        os.replace(bb + ".tmp", bb)
                elif m == "progdp":
                    if not glob.glob(os.path.join(bb, "*")):
                        raise RuntimeError("progdp needs the gcm backbones; run gcm first")
                    mergers.gcm(a, b, out, mwork, 1, backbones=bb, trace="progdp")
                rec["seconds"] = round(time.time() - start, 1)
                rec["constraintsKept"] = mergers.check_constraints(a, b, out)
                est = fasta.upper(fasta.read(out))
                rec.update(score.fastsp(os.path.join(work, "true.fasta"), out))
                rec.update(cross_scores(true, est, A))
            except Exception as exc:  # record and continue with the other mergers
                rec["error"] = repr(exc)[:300]
                rec["seconds"] = round(time.time() - start, 1)
            records.append(rec)
            print(json.dumps(rec), flush=True)
            with open(results, "a") as f:  # one short line per append: safe across worker processes
                f.write(json.dumps(rec) + "\n")
    return records


def main():
    p = argparse.ArgumentParser()
    p.add_argument("results")
    p.add_argument("replicates", nargs="+", help="DATASET:REP_DIR")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--conditions", default="oracle,fftnsi")
    p.add_argument("--mergers", default=",".join(ALL_MERGERS))
    args = p.parse_args()
    done = set()
    if os.path.exists(args.results):
        with open(args.results) as f:
            done = {(r["dataset"], r["condition"], r["merger"]) for r in map(json.loads, f) if "avgErr" in r}
    conditions = args.conditions.split(",")
    names = args.mergers.split(",")
    with ProcessPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_replicate, args.results, *r.split(":", 1), conditions, names, done) for r in args.replicates]
        for fut in futures:
            try:
                fut.result()
            except Exception as exc:
                print("replicate failed:", repr(exc), file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
