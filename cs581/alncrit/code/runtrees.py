"""FastTree-2 (GTR+G, nucleotides) on alignments, scored against the true tree.

    python runtrees.py JOBS.tsv OUT.jsonl [--workers 2] [--method ft|iqfast]

JOBS.tsv lines: KEY  ALIGNMENT  TRUE_TREE  TREE_OUT   (KEY is e.g. 1000M2/R0/gcm)
Rows already in OUT.jsonl (by KEY+method) are skipped, so the script is restartable.
"""
import argparse
import json
import os
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import treeerr

IQTREE = "/opt/mm/root/envs/bio/bin/iqtree3"


def run(job, method):
    key, aln, true_tree, out = job
    os.makedirs(os.path.dirname(out), exist_ok=True)
    start = time.time()
    if not os.path.exists(out):
        if method == "ft":
            with open(out + ".part", "w") as f:
                subprocess.run(["FastTree", "-nt", "-gtr", "-gamma", "-quiet", "-nosupport", aln],
                               stdout=f, stderr=subprocess.DEVNULL, check=True)
            os.rename(out + ".part", out)
        else:  # IQ-TREE --fast, GTR+G4, 1 thread
            with tempfile.TemporaryDirectory(dir="/opt/runs") as tmp:
                subprocess.run([IQTREE, "-s", aln, "-m", "GTR+G4", "--fast", "-T", "1", "-pre", tmp + "/x",
                                "--quiet", "-redo"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                shutil.copy(tmp + "/x.treefile", out)
    secs = time.time() - start
    row = {"key": key, "method": method, "seconds": round(secs, 1), "aln": aln}
    row.update(treeerr.error(true_tree, out))
    return row


def main():
    p = argparse.ArgumentParser()
    p.add_argument("jobs")
    p.add_argument("out")
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--method", default="ft")
    a = p.parse_args()
    done = set()
    if os.path.exists(a.out):
        done = {(r["key"], r["method"]) for r in map(json.loads, open(a.out))}
    jobs = [tuple(l.split()) for l in open(a.jobs) if l.strip()]
    jobs = [j for j in jobs if (j[0], a.method) not in done]
    print(len(jobs), "jobs", flush=True)
    with ProcessPoolExecutor(a.workers) as pool:
        futs = {pool.submit(run, j, a.method): j for j in jobs}
        for f in as_completed(futs):
            try:
                row = f.result()
            except Exception as e:
                print("FAILED", futs[f][0], repr(e)[:200], flush=True)
                continue
            with open(a.out, "a") as o:
                o.write(json.dumps(row) + "\n")
            print(row["key"], row["method"], row["seconds"], round(row["fn_rate"], 4), flush=True)


if __name__ == "__main__":
    main()
