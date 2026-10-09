"""Run ML tree estimators on cached alignments and score them against the true tree.

    python runtrees.py OUT.jsonl --datasets 1000M2 1000M3 --reps 0 1 2 \
        --alns true_align gcm pasta_align --methods fasttree iqtree raxmlng [--workers 4]

Each job uses one thread; `--workers` jobs run concurrently. Finished (dataset, rep, aln,
method) rows already present in OUT.jsonl are skipped, so runs are resumable.
Trees are kept in $MLDATA/<ds>/R<rep>/trees/<aln>.<method>.tre.
"""
import argparse
import json
import re
import traceback
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import treeerr  # noqa: E402

MLDATA = os.environ.get("MLDATA", "/opt/data/mlcache")
BIO = "/opt/mm/root/envs/bio/bin"
FASTTREE = shutil.which("FastTree") or BIO + "/fasttree"
IQTREE = BIO + "/iqtree3"
RAXMLNG = shutil.which("raxml-ng") or BIO + "/raxml-ng"


def read_fasta(path):
    names, seqs, cur = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if names:
                    seqs.append("".join(cur))
                names.append(line[1:].split()[0])
                cur = []
            elif line:
                cur.append(line)
    seqs.append("".join(cur))
    return names, seqs


def write_fasta(path, names, seqs):
    with open(path, "w") as f:
        for n, s in zip(names, seqs):
            f.write(">%s\n%s\n" % (n, s))


def clean_alignment(src, dst):
    """Upper-case, map '.'/'?' to '-', drop all-gap columns (some tools reject them)."""
    names, seqs = read_fasta(src)
    seqs = [s.upper().replace(".", "-").replace("?", "-") for s in seqs]
    keep = [i for i in range(len(seqs[0])) if any(s[i] != "-" for s in seqs)]
    write_fasta(dst, names, ["".join(s[i] for i in keep) for s in seqs])
    return len(keep)


def run(cmd, log, cwd=None, stdout=None):
    t = time.time()
    with open(log, "w") as lf:
        subprocess.run(cmd, cwd=cwd, stdout=stdout or lf, stderr=lf, check=True)
    return time.time() - t


def estimate(method, aln, out_tree, work, extra=None):
    """Run one estimator; returns (seconds, log-likelihood reported by the tool or None)."""
    extra = extra or []
    if method.startswith("fasttree"):
        opts = {"fasttree": ["-gtr", "-gamma"], "fasttree_fastest": ["-gtr", "-gamma", "-fastest"],
                "fasttree_pasta": ["-gtr", "-gamma", "-fastest"]}[method]
        with open(out_tree, "w") as o:
            sec = run([FASTTREE, "-nt"] + opts + extra + [aln], os.path.join(work, "ft.log"), stdout=o)
        lnl = None
        for line in open(os.path.join(work, "ft.log")):
            m = re.match(r"Gamma\(20\) LogLk = (\S+)", line)
            if m:
                lnl = float(m.group(1))
        return sec, lnl
    if method.startswith("iqtree"):
        opts = {"iqtree": [], "iqtree_fast": ["--fast"]}[method]
        sec = run([IQTREE, "-s", aln, "-m", "GTR+G", "-T", "1", "--seed", "1", "--prefix",
                   os.path.join(work, "iq"), "-redo", "--quiet"] + opts + extra, os.path.join(work, "iq.out"))
        shutil.copy(os.path.join(work, "iq.treefile"), out_tree)
        lnl = None
        for line in open(os.path.join(work, "iq.log")):
            if line.startswith("BEST SCORE FOUND"):
                lnl = float(line.split(":")[1])
        return sec, lnl
    if method.startswith("raxmlng"):
        start = {"raxmlng": ["--tree", "pars{1}"], "raxmlng_rand1": ["--tree", "rand{1}"],
                 "raxmlng_default": []}[method]
        sec = run([RAXMLNG, "--search", "--msa", aln, "--model", "GTR+G", "--threads", "1", "--seed", "1",
                   "--prefix", os.path.join(work, "rx"), "--redo"] + start + extra, os.path.join(work, "rx.out"))
        shutil.copy(os.path.join(work, "rx.raxml.bestTree"), out_tree)
        lnl = None
        for line in open(os.path.join(work, "rx.raxml.log")):
            if line.startswith("Final LogLikelihood:"):
                lnl = float(line.split(":")[1])
        return sec, lnl
    raise ValueError(method)


def job(ds, rep, aln, method, out_path, lock):
    d = os.path.join(MLDATA, ds, "R%s" % rep)
    os.makedirs(os.path.join(d, "trees"), exist_ok=True)
    clean = os.path.join(d, aln + ".clean.fasta")
    if not os.path.exists(clean):
        clean_alignment(os.path.join(d, aln + ".fasta"), clean + ".tmp")
        os.replace(clean + ".tmp", clean)
    tree = os.path.join(d, "trees", "%s.%s.tre" % (aln, method))
    work = tempfile.mkdtemp(prefix="ml_%s_%s_%s_%s_" % (ds, rep, aln, method))
    try:
        sec, lnl = estimate(method, clean, tree, work)
    except subprocess.CalledProcessError as e:
        print("FAILED", ds, rep, aln, method, e, "logs in", work, flush=True)
        return
    err = treeerr.error(os.path.join(d, "true_tree.tre"), tree)
    row = {"dataset": ds, "rep": rep, "aln": aln, "method": method, "seconds": round(sec, 1),
           "lnl_tool": lnl, **{k: err[k] for k in ("fn_rate", "fp_rate", "rf_rate")}, "tree": tree}
    with lock:
        with open(out_path, "a") as f:
            f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
    shutil.rmtree(work, ignore_errors=True)


def main():
    import threading
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--datasets", nargs="+", required=True)
    ap.add_argument("--reps", nargs="+", default=["0", "1", "2"])
    ap.add_argument("--alns", nargs="+", default=["true_align", "gcm", "pasta_align"])
    ap.add_argument("--methods", nargs="+", default=["fasttree", "iqtree", "raxmlng"])
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    done = set()
    if os.path.exists(a.out):
        for line in open(a.out):
            r = json.loads(line)
            done.add((r["dataset"], str(r["rep"]), r["aln"], r["method"]))
    jobs = [(ds, rep, aln, m) for m in a.methods for ds in a.datasets for rep in a.reps for aln in a.alns
            if (ds, rep, aln, m) not in done]
    print(len(jobs), "jobs", flush=True)
    lock = threading.Lock()
    with ThreadPoolExecutor(a.workers) as ex:
        futs = [ex.submit(job, *j, a.out, lock) for j in jobs]
        for f in futs:
            try:
                f.result()
            except Exception:
                traceback.print_exc()


if __name__ == "__main__":
    main()
