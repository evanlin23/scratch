"""MAGUS with different subset (base) aligners: paired, merge-only.

    PYTHONPATH=cs581/code python3 cs581/basemeth/code/basemeth.py prep  JOBS WORK
    PYTHONPATH=cs581/code python3 cs581/basemeth/code/basemeth.py align JOBS WORK --methods a,b [--jobs 4]
    PYTHONPATH=cs581/code python3 cs581/basemeth/code/basemeth.py merge JOBS WORK OUT.jsonl --methods a,b

JOBS lines: NAME TRUE_ALIGNMENT NUM_SUBSETS [UNALIGNED | cache:RUN]. With UNALIGNED (HomFam) MAGUS aligns
that file and every score is on the estimate restricted to the reference (seed) sequences. With cache:RUN
the subsets and backbones come from cs581/experiments/runs/RUN/inputs.tar.xz (an earlier MAGUS run with the
same flags) instead of a new MAGUS run.

prep   one MAGUS run per dataset (paper flags, 4 threads, timed); keeps WORK/NAME/inputs/{subalignments,backbones}.
align  every subset re-aligned from its unaligned sequences by every method, one thread each, --jobs subsets in
       parallel; per-subset wall and CPU seconds and the subset's own FastSP score -> WORK/NAME/align.jsonl.
merge  the GCM merge only (gcmx.run_magus -s SUBSETS -b BACKBONES, 4 threads) on each method's subsets with
       MAGUS's own backbones; final FastSP (SPFN, SPFP, TC) -> OUT.jsonl. Restartable: finished rows are skipped.
"""

import argparse
import concurrent.futures
import json
import os
import resource
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time

from gcmx import fasta, score
from gcmx.e2e_bench import CODE, REPO, magus_flags, seq_type

MAGUS_MAFFT = os.path.join(CODE, "MAGUS", "magus", "tools", "mafft", "mafft")
BIO = "/opt/mm/root/envs/bio/bin"
ALN2 = "/opt/mm/root/envs/aln2/bin"
MERGE_FLAGS = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
               "-f", "4"]

# one thread each; {in} unaligned fasta, {out} output fasta, {outbase} prefix (Prank)
METHODS = {
    "linsi": [MAGUS_MAFFT, "--localpair", "--maxiterate", "1000", "--ep", "0.123", "--quiet", "--thread", "1",
              "--anysymbol", "{in}"],  # MAGUS's own subset command
    "ginsi": [MAGUS_MAFFT, "--globalpair", "--maxiterate", "1000", "--quiet", "--thread", "1", "--anysymbol", "{in}"],
    "muscle5": [ALN2 + "/muscle", "-align", "{in}", "-output", "{out}", "-threads", "1"],
    "famsa": [BIO + "/famsa", "-t", "1", "{in}", "{out}"],
    "probcons": [ALN2 + "/probcons", "{in}"],
    "clustalo": ["clustalo", "-i", "{in}", "-o", "{out}", "--threads", "1", "--force"],
    "kalign": [ALN2 + "/kalign", "-i", "{in}", "-o", "{out}", "-n", "1"],
    "prank": [ALN2 + "/prank", "-d={in}", "-o={outbase}", "-f=fasta", "-quiet"],
    "tcoffee": [ALN2 + "/t_coffee", "{in}", "-output", "fasta_aln", "-outfile", "{out}", "-quiet", "-n_core", "1"],
}


def jobs(path):
    out = []
    for line in open(path):
        f = line.split()
        if not f:
            continue
        true = f[1] if os.path.isabs(f[1]) else os.path.join(REPO, f[1])
        extra = f[3] if len(f) > 3 else None
        out.append({"name": f[0], "true": true, "k": f[2],
                    "unaligned": extra if extra and not extra.startswith("cache:") else None,
                    "cache": extra[6:] if extra and extra.startswith("cache:") else None})
    return out


def timed(cmd, log, cwd=CODE, stdout=None):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.time()
    with open(log, "w") as f:
        proc = subprocess.run(cmd, cwd=cwd, stdout=stdout or f, stderr=f if stdout else subprocess.STDOUT)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    if proc.returncode != 0:
        raise RuntimeError("failed ({}): {}".format(proc.returncode, " ".join(cmd)))
    return round(time.time() - start, 2), round(cpu, 2)


def accuracy(ref, est_path):
    """FastSP of est against ref, restricted to ref's sequences when est has more (HomFam)."""
    est = fasta.upper(fasta.read(est_path))
    with tempfile.TemporaryDirectory() as tmp:
        r, e = os.path.join(tmp, "r.fa"), os.path.join(tmp, "e.fa")
        fasta.write(ref, r)
        fasta.write(fasta.restrict(est, list(ref)) if len(est) != len(ref) else est, e)
        s = score.fastsp(r, e)
    return {k: s[k] for k in ("SPFN", "SPFP", "avgErr", "TC", "shared", "refHom", "estHom", "LenEst", "LenRef")}


def prep(args):
    for j in jobs(args.jobs):
        w = os.path.join(args.work, j["name"])
        inputs = os.path.join(w, "inputs")
        if os.path.exists(os.path.join(w, "prep.json")):
            continue
        os.makedirs(w, exist_ok=True)
        ref = fasta.upper(fasta.read(j["true"]))
        true = os.path.join(w, "true.fasta")
        fasta.write(ref, true)
        info = {"dataset": j["name"], "type": seq_type(ref), "nseq": None}
        if j["cache"]:
            shutil.rmtree(inputs, ignore_errors=True)
            with tarfile.open(os.path.join(REPO, "cs581", "experiments", "runs", j["cache"], "inputs.tar.xz")) as t:
                t.extractall(w)
            src = os.path.join(REPO, "cs581", "experiments", "runs", j["cache"], "prep.json")
            p = json.load(open(src))
            info.update({"magus_wall": p.get("seconds"), "magus_avgErr": p.get("avgErr"), "source": "cache:" + j["cache"]})
        else:
            unaligned = j["unaligned"] or os.path.join(w, "unaligned.fasta")
            if not j["unaligned"]:
                fasta.write(fasta.ungap(ref), unaligned)
            out = os.path.join(w, "magus.fasta")
            if os.path.exists(out):
                os.remove(out)
            shutil.rmtree(os.path.join(w, "magus"), ignore_errors=True)
            wall, cpu = timed([sys.executable, "-m", "gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", "4", "-d",
                               os.path.join(w, "magus"), "-i", unaligned, "-o", out] + magus_flags(j["k"]),
                              os.path.join(w, "magus.log"))
            shutil.rmtree(inputs, ignore_errors=True)
            shutil.copytree(os.path.join(w, "magus", "subalignments"), os.path.join(inputs, "subalignments"))
            os.makedirs(os.path.join(inputs, "backbones"))
            for f in os.listdir(os.path.join(w, "magus", "graph")):
                if f.startswith("backbone_") and f.endswith("_mafft.txt"):
                    shutil.copy(os.path.join(w, "magus", "graph", f), os.path.join(inputs, "backbones"))
            shutil.rmtree(os.path.join(w, "magus"), ignore_errors=True)
            a = accuracy(ref, out)
            info.update({"magus_wall": wall, "magus_cpu": cpu, "magus_avgErr": a["avgErr"], "magus": a,
                         "source": "magus"})
        info["nseq"] = sum(len(fasta.read(os.path.join(inputs, "subalignments", f)))
                           for f in os.listdir(os.path.join(inputs, "subalignments")))
        json.dump(info, open(os.path.join(w, "prep.json"), "w"))
        print(json.dumps(info), flush=True)


def run_method(method, unaligned, out_path, tmp):
    inp = os.path.join(tmp, "in.fa")
    shutil.copy(unaligned, inp)
    outbase = os.path.join(tmp, "out")
    out = outbase + ".fa"
    argv = METHODS[method]
    cmd = [a.format(**{"in": inp, "out": out, "outbase": outbase}) for a in argv]
    log = os.path.join(tmp, "log.txt")
    if any("{out}" in a or "{outbase}" in a for a in argv):
        wall, cpu = timed(cmd, log, cwd=tmp)
    else:
        with open(out, "w") as o:
            wall, cpu = timed(cmd, log, cwd=tmp, stdout=o)
    if method == "prank":
        out = outbase + ".best.fas"
    aln = fasta.upper(fasta.read(out))
    assert set(aln) == set(fasta.read(unaligned)), "missing sequences"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fasta.write(aln, out_path)
    return wall, cpu


def align_one(w, method, filename, ref):
    sub = fasta.read(os.path.join(w, "inputs", "subalignments", filename))
    out_path = os.path.join(w, "sub_" + method, filename)
    row = {"method": method, "subset": filename, "nseq": len(sub)}
    with tempfile.TemporaryDirectory() as tmp:
        unaligned = os.path.join(tmp, "unaligned.fa")
        fasta.write(fasta.upper(fasta.ungap(sub)), unaligned)
        try:
            row["wall"], row["cpu"] = run_method(method, unaligned, out_path, tmp)
        except Exception as e:
            row["error"] = repr(e)[:300]
            return row
    seeds = [n for n in sub if n in ref]
    if len(seeds) >= 2:
        r = {n: ref[n] for n in seeds}
        r = fasta.restrict(r, seeds)
        with tempfile.TemporaryDirectory() as tmp:
            rp, ep = os.path.join(tmp, "r.fa"), os.path.join(tmp, "e.fa")
            fasta.write(r, rp)
            est = fasta.upper(fasta.read(out_path))
            fasta.write(fasta.restrict(est, seeds) if len(seeds) != len(est) else est, ep)
            try:
                s = score.fastsp(rp, ep)
                row.update({k: s[k] for k in ("SPFN", "SPFP", "shared", "refHom", "estHom")})
            except Exception as e:  # e.g. a 2-seed subset with no homologies
                row["score_error"] = repr(e)[:100]
    return row


def align(args):
    methods = args.methods.split(",")
    tasks = []
    for j in jobs(args.jobs):
        w = os.path.join(args.work, j["name"])
        if not os.path.exists(os.path.join(w, "prep.json")):
            continue
        log = os.path.join(w, "align.jsonl")
        done = set()
        if os.path.exists(log):
            done = {(r["method"], r["subset"]) for r in map(json.loads, open(log)) if "error" not in r}
        ref = fasta.upper(fasta.read(os.path.join(w, "true.fasta")))
        for m in methods:
            for f in sorted(os.listdir(os.path.join(w, "inputs", "subalignments"))):
                if (m, f) not in done:
                    os.makedirs(os.path.join(w, "sub_" + m), exist_ok=True)
                    tasks.append((w, m, f, ref))
    print("align tasks:", len(tasks), flush=True)
    with concurrent.futures.ThreadPoolExecutor(args.jobs_parallel) as pool:
        futs = {pool.submit(align_one, *t): t for t in tasks}
        for fut in concurrent.futures.as_completed(futs):
            row = fut.result()
            w = futs[fut][0]
            with open(os.path.join(w, "align.jsonl"), "a") as f:
                f.write(json.dumps(row) + "\n")
            if "error" in row:
                print(os.path.basename(w), json.dumps(row), flush=True)


def merge(args):
    methods = args.methods.split(",")
    done = set()
    if os.path.exists(args.out):
        done = {(r["dataset"], r["method"]) for r in map(json.loads, open(args.out))}
    for j in jobs(args.jobs):
        w = os.path.join(args.work, j["name"])
        if not os.path.exists(os.path.join(w, "prep.json")):
            continue
        ref = fasta.upper(fasta.read(os.path.join(w, "true.fasta")))
        rows = {}
        log = os.path.join(w, "align.jsonl")
        for r in (map(json.loads, open(log)) if os.path.exists(log) else []):
            rows[(r["method"], r["subset"])] = r
        nsub = len(os.listdir(os.path.join(w, "inputs", "subalignments")))
        for m in methods:
            if (j["name"], m) in done:
                continue
            sub = [r for (mm, _), r in rows.items() if mm == m]
            if m == "orig":  # control: MAGUS's own subset alignments, unchanged
                sub = [{"wall": 0, "cpu": 0}] * nsub
            if len(sub) < nsub or any("error" in r for r in sub):
                continue
            d = os.path.join(w, "merge_" + m)
            shutil.rmtree(d, ignore_errors=True)
            out = os.path.join(w, "merge_" + m + ".fasta")
            if os.path.exists(out):
                os.remove(out)
            wall, cpu = timed([sys.executable, "-m", "gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", "4",
                               "-d", d, "-s", os.path.join(w, "inputs", "subalignments") if m == "orig" else os.path.join(w, "sub_" + m), "-b", os.path.join(w, "inputs", "backbones"),
                               "-o", out] + MERGE_FLAGS, os.path.join(w, "merge_" + m + ".log"))
            shutil.rmtree(d, ignore_errors=True)
            row = {"dataset": j["name"], "method": m, "merge_wall": wall, "merge_cpu": cpu,
                   "sub_wall": round(sum(r["wall"] for r in sub), 1), "sub_cpu": round(sum(r["cpu"] for r in sub), 1),
                   "sub_max_wall": max(r["wall"] for r in sub), **accuracy(ref, out)}
            scored = [r for r in sub if "shared" in r]
            if scored:
                sh, rh, eh = (sum(r[k] for r in scored) for k in ("shared", "refHom", "estHom"))
                row.update({"sub_SPFN": 1 - sh / rh if rh else None, "sub_SPFP": 1 - sh / eh if eh else None,
                            "sub_scored": len(scored)})
            with open(args.out, "a") as f:
                f.write(json.dumps(row) + "\n")
            print(json.dumps(row), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=("prep", "align", "merge"))
    p.add_argument("jobs")
    p.add_argument("work")
    p.add_argument("out", nargs="?")
    p.add_argument("--methods", default="linsi")
    p.add_argument("--jobs-parallel", type=int, default=4)
    args = p.parse_args()
    {"prep": prep, "align": align, "merge": merge}[args.cmd](args)


if __name__ == "__main__":
    main()
