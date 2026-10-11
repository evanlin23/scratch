# AI-assisted (Claude), exploration code for CS581 project
"""Regenerate every estimated alignment that has a tree row (data/trees.jsonl) from the bank replicates.

    python3 regen.py [--lanes 4] [--only METHOD,...]

Each job runs single-threaded (--threads 1 for vote.py; gg.py with bbe.THREADS = 1) so 4 lanes fill 4 cores.
Outputs: /opt/work/treecrit/reps/<key>/vote/<variant>/out.fasta (vote.py) or variants/<...>/out.fasta (gg.py).
Restartable: jobs whose output exists are skipped. data/aln_paths.jsonl maps (key, method) -> path.
"""
import argparse
import json
import os
import subprocess
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
CS = os.path.abspath(os.path.join(HERE, "..", ".."))
W = "/opt/work/treecrit/reps"
VOTE = {"magus": "magus", "es4": "es4", "vote_hard": "hard+mask", "vote_hard_mask": "hard+mask",
        "vote_hard-bb": "hard-bb", "vote_soft": "soft", "vote_soft-bb": "soft-bb", "vote_soft2": "soft2",
        "vote_soft4": "soft4"}
GG = {"recipe": "wsoft0.03:linsi&fftns2#es4", "es3": "linsi#es3", "hardf": "linsi&fftns2-op3"}
# h8 (SIMHIGH_R17) ran no vote.py es4; its es4 tree is on gg.py's linsi#es4 alignment
GG_OVERRIDE = {("SIMHIGH_R17_gv", "es4"): "linsi#es4"}
PRIO = ["magus", "es4", "vote_hard-bb", "recipe", "hardf", "es3", "vote_hard", "vote_soft", "vote_soft-bb",
        "vote_soft2", "vote_soft4"]


def gg_dir(v):
    safe = v.replace("&", "_i_").replace("|", "_m_").replace("+", "_p_").replace("~", "_k_").replace("^", "_l_")
    return safe.replace(":", "_c_").replace("#", "_es_")  # = gg.py's directory name (bbe.safe + gg)


def path_of(key, m):
    rep = os.path.join(W, key)
    if m == "true":
        return os.path.join(rep, "true.fasta")
    if (key, m) in GG_OVERRIDE:
        return os.path.join(rep, "variants", gg_dir(GG_OVERRIDE[(key, m)]), "out.fasta")
    if m in VOTE:
        v = VOTE[m]
        return os.path.join(rep, "vote", v, "out.masked.fasta" if m == "vote_hard_mask" else "out.fasta")
    if m in GG:
        return os.path.join(rep, "variants", gg_dir(GG[m]), "out.fasta")
    return None


def job(arg):
    key, m = arg
    rep = os.path.join(W, key)
    out = path_of(key, m)
    if os.path.exists(out) and os.path.exists(os.path.join(os.path.dirname(out), "done")):
        return key, m, "skip"
    ggv = GG_OVERRIDE.get((key, m), GG.get(m))
    if m in VOTE and (key, m) not in GG_OVERRIDE:
        cmd = [sys.executable, os.path.join(CS, "gcmvote", "code", "vote.py"), rep, VOTE[m], "--threads", "1"]
    else:
        lock = os.path.join(rep, "variants", gg_dir(ggv) + ".lock")
        if os.path.exists(lock):
            os.remove(lock)  # stale lock from an interrupted run
        cmd = [sys.executable, "-c", "import sys; sys.path.insert(0, %r); import bbe; bbe.THREADS = 1; "
               "sys.argv = ['gg.py', 'run', %r, %r]; sys.path.insert(0, %r); import runpy; "
               "runpy.run_path(%r, run_name='__main__')" % (os.path.join(CS, "bbevidence", "code"), rep, ggv,
                                                           os.path.join(CS, "gcmgen", "code"),
                                                           os.path.join(CS, "gcmgen", "code", "gg.py"))]
    log = open(os.path.join(rep, "regen_{}.log".format(m)), "w")
    r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=os.path.join(CS, "code"))
    if r.returncode == 0 and os.path.exists(out):
        open(os.path.join(os.path.dirname(out), "done"), "w").close()
        return key, m, "ok"
    return key, m, "FAILED rc={}".format(r.returncode)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lanes", type=int, default=4)
    ap.add_argument("--only")
    ap.add_argument("--paths-only", action="store_true")
    a = ap.parse_args()
    if a.paths_only:
        return write_paths()
    trees = [json.loads(l) for l in open(os.path.join(HERE, "..", "data", "trees.jsonl"))]
    need = sorted({(r["key"], r["method"]) for r in trees if r["method"] in VOTE or r["method"] in GG},
                  key=lambda x: (PRIO.index(x[0 + 1]) if x[1] in PRIO else 99, x[0]))
    if a.only:
        need = [x for x in need if x[1] in a.only.split(",")]
    # hard+mask serves both vote_hard and vote_hard_mask: run it once
    jobs, seen = [], set()
    for key, m in need:
        k = (key, VOTE.get(m, GG.get(m)))
        if k in seen:
            continue
        seen.add(k)
        jobs.append((key, m))
    print(len(jobs), "jobs", flush=True)
    with Pool(a.lanes) as p:
        for res in p.imap_unordered(job, jobs):
            print(*res, flush=True)
    write_paths()


def write_paths():
    trees = [json.loads(l) for l in open(os.path.join(HERE, "..", "data", "trees.jsonl"))]
    with open(os.path.join(HERE, "..", "data", "aln_paths.jsonl"), "w") as f:
        for r in trees:
            p = path_of(r["key"], r["method"])
            if p and os.path.exists(p) and os.path.exists(os.path.join(os.path.dirname(p), "done")):
                f.write(json.dumps({"key": r["key"], "method": r["method"], "path": p}) + "\n")


if __name__ == "__main__":
    main()
