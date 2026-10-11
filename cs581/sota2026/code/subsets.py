"""Subset-level test: which aligner beats MAFFT L-INS-i on MAGUS-sized subsets?

    python3 subsets.py plan PLAN.jsonl            # build the subset list (deterministic)
    python3 subsets.py run PLAN.jsonl OUT.jsonl --tools a,b --workers 4 [--kinds magus40,rand40] [--max-n 60]

Subset kinds, per replicate:
  magus40  2 of MAGUS's own 25 decomposition subsets (centroid-edge split of its guide tree;
           ~40 sequences on 1000-sequence data), read from the cached MAGUS run in
           cs581/experiments/runs/<rep>/inputs.tar.xz. Baseline "cached-linsi" = the L-INS-i
           alignment MAGUS itself produced for that subset.
  bb200    MAGUS's first backbone (8 random sequences per subset, 200 total), same cache.
  rand40 / rand200   uniformly random sequences from the full reference (not tree-local:
           harder than what MAGUS gives its base method).
Every subset is scored with FastSP against the reference alignment restricted to it.
Jobs run single-threaded, `workers` at a time (accuracy is the target here; times are 1-core).
"""

import argparse
import concurrent.futures as cf
import json
import os
import random
import shutil
import subprocess
import sys
import tarfile
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "cs581", "code"))
sys.path.insert(0, HERE)
from gcmx import fasta, score  # noqa: E402
import tools  # noqa: E402

RUNS = os.path.join(REPO, "cs581", "experiments", "runs")
WORK = "/opt/runs/sota/subsets"
D = "/opt/data/Datasets"


def reference(rep):
    name = rep.rsplit("_", 1)[0]
    if name.startswith("1000"):
        return "{}/ROSE/{}/R0/rose.aln.true.fasta".format(D, name)
    if name == "RNASim":
        return D + "/RNASim/1000/R0/true_align.txt"
    if name.startswith("BBA"):
        return os.path.join(REPO, "cs581/data/balibase_clean/RV100_{}.fasta".format(name))
    if name.startswith("16S"):
        return "{}/Gutell/{}/R0/true_align_clean.txt".format(D, name)
    raise ValueError(rep)


def cached(rep):
    """Extract a replicate's cached MAGUS inputs once; returns the inputs/ dir or None."""
    tar = os.path.join(RUNS, rep, "inputs.tar.xz")
    if not os.path.exists(tar):
        return None
    dest = os.path.join(WORK, "cache", rep)
    if not os.path.isdir(os.path.join(dest, "inputs")):
        os.makedirs(dest, exist_ok=True)
        with tarfile.open(tar) as t:
            t.extractall(dest, filter="data")
    return os.path.join(dest, "inputs")


def plan(out):
    reps = ["1000L1_R0", "1000L2_R0", "1000L3_R0", "1000M1_R0", "1000M2_R0", "1000M3_R0", "1000M4_R0", "1000S1_R0",
            "1000S2_R0", "1000S3_R0", "RNASim_R0", "16S.M_R0", "BBA0067_R0", "BBA0101_R0", "BBA0154_R0", "BBA0190_R0"]
    rows = []
    for rep in reps:
        rng = random.Random(rep)
        names = list(fasta.read(reference(rep)))
        inp = cached(rep)
        if inp:
            subs = sorted(os.listdir(os.path.join(inp, "subalignments")), key=lambda s: int(s.split("_")[-1][:-4]))
            for s in rng.sample(subs, 2):
                p = os.path.join(inp, "subalignments", s)
                rows.append({"rep": rep, "kind": "magus40", "id": s[:-4], "taxa": list(fasta.read(p)), "cached": p})
            p = os.path.join(inp, "backbones", "backbone_1_mafft.txt")
            rows.append({"rep": rep, "kind": "bb200", "id": "backbone_1", "taxa": list(fasta.read(p)), "cached": p})
        rows.append({"rep": rep, "kind": "rand40", "id": "rand40_0", "taxa": rng.sample(names, 40)})
        rows.append({"rep": rep, "kind": "rand200", "id": "rand200_0", "taxa": rng.sample(names, min(200, len(names)))})
    with open(out, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    print(len(rows), "subsets")


lock = threading.Lock()


def one(task, tool, out_path, timeout):
    w = os.path.join(WORK, "jobs", task["rep"], task["id"])
    os.makedirs(w, exist_ok=True)
    true, unal = os.path.join(w, "true.fasta"), os.path.join(w, "unaligned.fasta")
    with lock:
        if not os.path.exists(unal):
            ref = fasta.upper(fasta.read(reference(task["rep"])))
            sub = fasta.restrict(ref, task["taxa"])
            fasta.write(sub, true)
            fasta.write(fasta.ungap(sub), unal)
    row = {"rep": task["rep"], "kind": task["kind"], "id": task["id"], "n": len(task["taxa"]), "tool": tool}
    est = os.path.join(w, tool + ".fasta")
    if tool == "cached-linsi":
        shutil.copy(task["cached"], est)
        row["status"] = "ok"
    else:
        tw = os.path.join(w, tool + "_work")
        os.makedirs(tw, exist_ok=True)
        row.update(tools.run(tool, unal, est, 1, tw, timeout, os.path.join(w, tool + ".log")))
        shutil.rmtree(tw, ignore_errors=True)
    if row["status"] == "ok":
        try:
            s = score.fastsp(true, est)
            row.update({k: s[k] for k in ("SPFN", "SPFP", "avgErr", "TC", "LenEst", "LenRef")})
        except Exception as e:
            row["status"] = "score failed: " + repr(e)[:200]
    with lock:
        with open(out_path, "a") as f:
            f.write(json.dumps(row) + "\n")
    print(json.dumps({k: row.get(k) for k in ("rep", "id", "tool", "status", "wall", "avgErr")}), flush=True)


def run(a):
    tasks = [json.loads(l) for l in open(a.plan)]
    kinds = set(a.kinds.split(","))
    done = set()
    if os.path.exists(a.out):
        done = {(r["rep"], r["id"], r["tool"]) for r in map(json.loads, open(a.out))}
    todo = []
    for t in tasks:
        if t["kind"] not in kinds or len(t["taxa"]) > a.max_n or (a.reps and t["rep"] not in a.reps.split(",")):
            continue
        for tool in a.tools.split(","):
            if (t["rep"], t["id"], tool) in done or (tool == "cached-linsi" and "cached" not in t):
                continue
            todo.append((t, tool))
    print(len(todo), "jobs", flush=True)
    with cf.ThreadPoolExecutor(a.workers) as ex:
        for f in [ex.submit(one, t, tool, a.out, a.timeout) for t, tool in todo]:
            f.result()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("cmd")
    p.add_argument("plan")
    p.add_argument("out", nargs="?")
    p.add_argument("--tools", default="")
    p.add_argument("--kinds", default="magus40,bb200,rand40,rand200")
    p.add_argument("--reps", default="")
    p.add_argument("--max-n", type=int, default=10 ** 9)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--timeout", type=float, default=1800)
    a = p.parse_args()
    plan(a.plan) if a.cmd == "plan" else run(a)
