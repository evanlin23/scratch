"""Run one method on one dataset made by make_long.py, score it, append JSON.

    python run.py DSDIR METHOD [--threads 4]

Common backbone (for the "add" methods, as UPP/WITCH would pick it):
sequences within 25% of the median length, with their TRUE alignment
(so only the adding step is tested), tree by FastTree; every other sequence
is a query.  Methods:

  magus        MAGUS on all sequences
  mafft        MAFFT --auto on all sequences
  upp          UPP (SEPP 4.4) -a backbone -t tree
  witch        WITCH (witch-msa) -b backbone -e tree
  emma         EMMA -b backbone -e tree
  mafft-add    mafft --add queries onto the backbone
  mafft-addlong mafft --addlong (MAFFT's mode for queries longer than the backbone)
  upp-trim / witch-trim   see trim.py (pilot fix)

Output: DSDIR/<method>/aln.fasta, DSDIR/<method>/score.json.
Scores: FastSP over all sequences, over the long sequences only and over the
others; for alignments with lower-case insertion letters also with -ml
(lower-case letters masked = treated as unaligned).
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))
from gcmx import fasta  # noqa: E402

ENV = "/opt/mm/root/envs/add/bin"
MAGUS = ["magus"]
FASTSP = "/opt/tools/FastSP/FastSP.jar"


def sh(cmd, log, **kw):
    with open(log, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True, **kw)


def fastsp(ref, est, mask_lower=False):
    cmd = ["java", "-Xmx4g", "-jar", FASTSP, "-r", ref, "-e", est] + (["-ml"] if mask_lower else [])
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    st = {}
    for line in (out.stdout + out.stderr).splitlines():
        t = line.split()
        if len(t) == 2 and t[0] in ("SPFN", "SPFP", "TC", "Compression"):
            st[t[0]] = float(t[1])
    st["err"] = (st["SPFN"] + st["SPFP"]) / 2
    return st


def score(ds, est_path):
    ref = fasta.read(os.path.join(ds, "true.fasta"))
    est = fasta.read(est_path)
    assert set(ref) == set(est), "taxa mismatch"
    long_ = [n for n in open(os.path.join(ds, "long.txt")).read().split() if n]
    rest = [n for n in ref if n not in set(long_)]
    # UPP/WITCH mark unaligned insertion letters in lower case; MAFFT writes all
    # lower case, so only mask when both cases occur
    text = "".join(est.values())
    has_lower = any(c.islower() for c in text) and any(c.isupper() for c in text)
    res = {}
    with tempfile.TemporaryDirectory() as tmp:
        for part, taxa in (("all", list(ref)), ("long", long_), ("rest", rest)):
            if len(taxa) < 2:
                continue
            r, e = os.path.join(tmp, "r.fa"), os.path.join(tmp, "e.fa")
            fasta.write(fasta.upper(fasta.restrict(ref, taxa)), r)
            sub = fasta.restrict(est, taxa)
            fasta.write(fasta.upper(sub), e)
            res[part] = fastsp(r, e)
            if has_lower:
                fasta.write(sub, e)
                res[part + "_ml"] = fastsp(r, e, mask_lower=True)
    return res


def backbone(ds, threads):
    """Backbone = sequences within 25% of the median length (true alignment) + FastTree tree."""
    bdir = os.path.join(ds, "backbone")
    done = os.path.join(bdir, "tree.nwk")
    os.makedirs(bdir, exist_ok=True)
    with open(os.path.join(bdir, ".lock"), "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)  # parallel jobs on one dataset share the backbone
        if not os.path.exists(done):
            make_backbone(ds, bdir, done)
    return bdir


def make_backbone(ds, bdir, done):
    seqs = fasta.read(os.path.join(ds, "unaligned.fasta"))
    med = statistics.median(len(s) for s in seqs.values())
    # queries = the designated (lengthened) sequences plus anything outside 25% of
    # the median, so a mult=0 control has exactly the same queries
    long_ = set(open(os.path.join(ds, "long.txt")).read().split())
    bb = {n: s for n, s in seqs.items() if 0.75 * med <= len(s) <= 1.25 * med and n not in long_}
    q = {n: s for n, s in seqs.items() if n not in bb}
    fasta.write(bb, os.path.join(bdir, "bb.unaln.fasta"))
    fasta.write(q, os.path.join(bdir, "queries.fasta"))
    # true alignment induced on the backbone sequences (isolates the adding step;
    # valid because no lengthened sequence is in the backbone)
    t0 = time.time()
    if os.environ.get("BACKBONE") == "magus":  # estimated backbone (validation runs)
        shutil.rmtree(os.path.join(bdir, "magus_wd"), ignore_errors=True)
        sh(MAGUS + ["-i", os.path.join(bdir, "bb.unaln.fasta"), "-o", os.path.join(bdir, "bb.aln.fasta"),
                    "-d", os.path.join(bdir, "magus_wd"), "-np", os.environ.get("THREADS", "4")],
           os.path.join(bdir, "log.txt"))
    else:
        ref = fasta.read(os.path.join(ds, "true.fasta"))
        fasta.write(fasta.restrict(ref, list(bb)), os.path.join(bdir, "bb.aln.fasta"))
    with open(os.path.join(bdir, "tree.tmp"), "w") as f:
        subprocess.run(["FastTree", "-nt", "-gtr", "-quiet", os.path.join(bdir, "bb.aln.fasta")],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)
    os.rename(os.path.join(bdir, "tree.tmp"), done)
    json.dump({"time": time.time() - t0, "n_backbone": len(bb), "n_query": len(q), "median": med},
              open(os.path.join(bdir, "info.json"), "w"))
    shutil.rmtree(os.path.join(bdir, "magus_wd"), ignore_errors=True)
    return bdir


def run_method(ds, method, threads, wd, log):
    """Returns path of the output alignment (all sequences)."""
    unal = os.path.join(ds, "unaligned.fasta")
    out = os.path.join(wd, "aln.fasta")
    if method == "magus":
        sh(MAGUS + ["-i", unal, "-o", out, "-d", os.path.join(wd, "wd"), "-np", str(threads)], log)
        return out
    if method == "mafft":
        with open(out, "w") as f:
            subprocess.run(["mafft", "--auto", "--thread", str(threads), unal], stdout=f,
                           stderr=open(log, "a"), check=True)
        return out
    bdir = backbone(ds, threads)
    bb, tree, q = (os.path.join(bdir, x) for x in ("bb.aln.fasta", "tree.nwk", "queries.fasta"))
    if method.endswith("-trim"):
        import trim
        q = trim.trim_queries(bb, q, wd, log)
        aln = add(method[:-len("-trim")], bb, tree, q, wd, out, threads, log)
        return trim.restore(aln, wd, os.path.join(wd, "restored.fasta"))
    return add(method, bb, tree, q, wd, out, threads, log)


def add(method, bb, tree, q, wd, out, threads, log):
    if method in ("mafft-add", "mafft-addlong"):
        flag = "--add" if method == "mafft-add" else "--addlong"
        with open(out, "w") as f:
            subprocess.run(["mafft", flag, q, "--thread", str(threads), bb], stdout=f,
                           stderr=open(log, "a"), check=True)
        return out
    if method == "upp":
        sh([ENV + "/python", ENV + "/run_upp.py", "-s", q, "-a", bb, "-t", tree, "-x", str(threads),
            "-d", wd + "/", "-o", "upp", "-m", "dna"], log)
        return os.path.join(wd, "upp_alignment.fasta")
    if method == "witch":
        sh([ENV + "/python", ENV + "/witch.py", "-b", bb, "-e", tree, "-q", q, "-t", str(threads),
            "-d", os.path.join(wd, "w"), "-o", "aln.fasta", "--molecule", "dna"], log)
        return os.path.join(wd, "w", "aln.fasta")
    if method == "emma":
        sh([ENV + "/python", "/opt/src/EMMA/emma.py", "-b", bb, "-e", tree, "-q", q, "-t", str(threads),
            "-d", os.path.join(wd, "e"), "-o", "aln.fasta", "--molecule", "dna"], log)
        return os.path.join(wd, "e", "aln.fasta")
    raise ValueError(method)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ds")
    ap.add_argument("method")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    wd = os.path.join(a.ds, a.method)
    res_path = os.path.join(wd, "score.json")
    if os.path.exists(res_path):
        print(open(res_path).read())
        return
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    log = os.path.join(wd, "log.txt")
    if a.method not in ("magus", "mafft"):
        backbone(a.ds, a.threads)  # not timed with the method; reported separately
    t0 = time.time()
    out = run_method(a.ds, a.method, a.threads, wd, log)
    el = time.time() - t0
    aln = fasta.read(out)
    # add-methods return backbone + queries: all sequences must be present
    res = {"ds": a.ds, "method": a.method, "time": el, "score": score(a.ds, out)}
    if a.method not in ("magus", "mafft"):
        res["backbone"] = json.load(open(os.path.join(a.ds, "backbone", "info.json")))
    shutil.copy(out, os.path.join(wd, "final.fasta")) if out != os.path.join(wd, "final.fasta") else None
    for d in ("wd", "w", "e"):
        shutil.rmtree(os.path.join(wd, d, "tree_decomp"), ignore_errors=True)
    json.dump(res, open(res_path, "w"))
    print(json.dumps(res))
    del aln


if __name__ == "__main__":
    main()
