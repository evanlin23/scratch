"""Aligners, FastTree, HMM add-back (UPP-style) and the three reference-free rogue detectors."""

import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GCMX_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "code"))
sys.path.insert(0, GCMX_ROOT)
from gcmx import fasta, score  # noqa: E402
from gcmx.e2e_bench import magus_flags, PASTA_PY, PASTA  # noqa: E402

TREESHRINK = "/opt/src/TreeShrink/run_treeshrink.py"
THREADS = "4"


def run(cmd, log, cwd=None, stdout=None):
    with open(log, "w") as lf:
        p = subprocess.run(cmd, cwd=cwd, stdout=stdout or lf, stderr=lf if stdout else subprocess.STDOUT)
    if p.returncode != 0:
        raise RuntimeError("failed: " + " ".join(cmd) + " (log " + log + ")")


def align(method, unaligned, work):
    """Align `unaligned` with method; returns (alignment dict, wall seconds)."""
    os.makedirs(work, exist_ok=True)
    out = os.path.join(work, method + ".fasta")
    log = os.path.join(work, method + ".log")
    t = time.time()
    if method == "mafft":
        with open(out, "w") as f:
            run(["mafft", "--auto", "--thread", THREADS, unaligned], log, stdout=f)
    elif method == "magus":
        d = os.path.join(work, "magus_tmp")
        shutil.rmtree(d, ignore_errors=True)
        run([sys.executable, "-m", "gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", THREADS, "-d", d,
             "-i", unaligned, "-o", out, "--graphbuildhmmextend", "false"] + magus_flags(25), log, cwd=GCMX_ROOT)
        shutil.rmtree(d, ignore_errors=True)
    elif method == "pasta":
        d, tmp = os.path.join(work, "pasta_out"), os.path.join(work, "pasta_tmp")
        for x in (d, tmp):
            shutil.rmtree(x, ignore_errors=True)
        os.makedirs(d)
        run([PASTA_PY, PASTA, "-i", unaligned, "-o", d, "-d", "dna", "--num-cpus", THREADS, "--iter-limit", "3",
             "--temporaries", tmp, "-j", "pastajob"], log)
        alns = [p for p in os.listdir(d) if p.startswith("pastajob.marker001.") and p.endswith(".aln") and "masked" not in p]
        shutil.copy(os.path.join(d, alns[0]), out)
        shutil.rmtree(d, ignore_errors=True)
        shutil.rmtree(tmp, ignore_errors=True)
    elif method == "linsi":
        with open(out, "w") as f:
            run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", THREADS, unaligned], log, stdout=f)
    else:
        raise ValueError(method)
    return fasta.upper(fasta.read(out)), round(time.time() - t, 1)


def fasttree(aln, path):
    """FastTree -nt -gtr on alignment dict; writes path (newick)."""
    fa = path + ".fa"
    fasta.write(aln, fa)
    with open(path, "w") as f:
        run(["FastTree", "-nt", "-gtr", "-quiet", "-nopr", "-nosupport", fa], path + ".log", stdout=f)
    os.remove(fa)
    return path


def drop_allgap(aln):
    names = list(aln)
    M = np.array([list(aln[n]) for n in names])
    keep = (M != "-").any(axis=0)
    return {n: "".join(M[i, keep]) for i, n in enumerate(names)}


def hmm_add(backbone, queries, work):
    """UPP-style: every backbone column is a match state (hmmbuild --hand); queries hmmaligned
    one at a time; query insertions get their own columns. Backbone homologies are unchanged."""
    os.makedirs(work, exist_ok=True)
    if not queries:
        return dict(backbone)
    names = list(backbone)
    L = len(backbone[names[0]])
    sto = os.path.join(work, "bb.sto")
    with open(sto, "w") as f:
        f.write("# STOCKHOLM 1.0\n")
        for n in names:
            f.write("%s %s\n" % (n, backbone[n].replace(".", "-")))
        f.write("#=GC RF %s\n//\n" % ("x" * L))
    hmm = os.path.join(work, "bb.hmm")
    run(["hmmbuild", "--hand", "--dna", "--cpu", THREADS, hmm, sto], os.path.join(work, "hmmbuild.log"))
    q = os.path.join(work, "q.fa")
    fasta.write(queries, q)
    a2m = os.path.join(work, "q.a2m")
    run(["hmmalign", "--dna", "--outformat", "A2M", "-o", a2m, hmm, q], os.path.join(work, "hmmalign.log"))
    qa = fasta.read(a2m)
    # per query: match-state chars (L of them) and insertion strings after each match state (index -1..L-1)
    parsed = {}
    for n, s in qa.items():
        match, ins, cur = [], {}, -1
        for ch in s:
            if ch.isupper() or ch == "-":
                match.append(ch)
                cur += 1
            else:
                ins.setdefault(cur, []).append(ch.upper())
        assert len(match) == L, (n, len(match), L)
        parsed[n] = (match, ins)
    width = {}
    for _, ins in parsed.values():
        for c, v in ins.items():
            width[c] = max(width.get(c, 0), len(v))
    out = {}
    for n in names:
        row = ["-" * width.get(-1, 0)]
        s = backbone[n]
        for c in range(L):
            row.append(s[c])
            row.append("-" * width.get(c, 0))
        out[n] = "".join(row)
    for n, (match, ins) in parsed.items():
        row = []
        for c in range(-1, L):
            if c >= 0:
                row.append(match[c])
            v = "".join(ins.get(c, []))
            row.append(v + "-" * (width.get(c, 0) - len(v)))
        out[n] = "".join(row)
    for n in queries:
        assert out[n].replace("-", "") == queries[n].upper(), n
    return out


def fastsp_restricted(true, est, taxa, work):
    os.makedirs(work, exist_ok=True)
    a, b = os.path.join(work, "ref.fa"), os.path.join(work, "est.fa")
    fasta.write(fasta.restrict(true, taxa), a)
    fasta.write(fasta.restrict(est, taxa), b)
    s = score.fastsp(a, b)
    return {k: s[k] for k in ("SPFN", "SPFP", "avgErr")}


# ---------------------------------------------------------------- detectors

def robust_z(x):
    x = np.asarray(x, float)
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    return (x - med) / (mad if mad > 0 else 1e-9)


def detect_treeshrink(tree_path, work):
    os.makedirs(work, exist_ok=True)
    out = os.path.join(work, "ts")
    shutil.rmtree(out, ignore_errors=True)
    run([sys.executable, "-W", "ignore", TREESHRINK, "-t", tree_path, "-o", out], os.path.join(work, "treeshrink.log"))
    # output.txt holds the removed leaves (output_summary.txt lists all candidates; not used)
    return set(open(os.path.join(out, "output.txt")).read().split())


def pdist_matrix(aln):
    names = list(aln)
    M = np.frombuffer("".join(aln[n] for n in names).encode(), dtype=np.uint8).reshape(len(names), -1)
    G = (M != ord("-")).astype(np.float32)
    same = np.zeros((len(names), len(names)), np.float32)
    for b in b"ACGT":
        X = (M == b).astype(np.float32)
        same += X @ X.T
    both = G @ G.T
    P = 1.0 - same / np.maximum(both, 1)
    P[both < 50] = 0.75
    return names, P


def detect_pdist(aln, nn=5, z=3.5):
    """Mean p-distance to the nn nearest neighbours; flag robust z > z."""
    names, P = pdist_matrix(aln)
    np.fill_diagonal(P, np.inf)
    score_ = np.sort(P, axis=1)[:, :nn].mean(axis=1)
    zz = robust_z(score_)
    return {n for n, v in zip(names, zz) if v > z}, dict(zip(names, score_.tolist()))


def detect_hmm(unaligned, work, nbb=100, z=3.5, seed=0):
    """UPP-style backbone: nbb random seqs with length within 25% of the median, MAFFT-L-INS-i,
    hmmbuild, hmmsearch every sequence; score = bit score per residue; flag robust z < -z."""
    import random
    os.makedirs(work, exist_ok=True)
    lens = {n: len(s) for n, s in unaligned.items()}
    med = float(np.median(list(lens.values())))
    pool = sorted(n for n, l in lens.items() if abs(l - med) <= 0.25 * med)
    bb = random.Random(seed).sample(pool, min(nbb, len(pool)))
    bbu = os.path.join(work, "bb_unaligned.fa")
    fasta.write({n: unaligned[n] for n in bb}, bbu)
    bba, _ = align("linsi", bbu, work)
    afa = os.path.join(work, "bb.afa")
    fasta.write(bba, afa)
    hmm = os.path.join(work, "bb.hmm")
    run(["hmmbuild", "--dna", "--cpu", THREADS, hmm, afa], os.path.join(work, "hmmbuild.log"))
    allq = os.path.join(work, "all.fa")
    fasta.write(unaligned, allq)
    tbl = os.path.join(work, "hits.tbl")
    run(["hmmsearch", "--max", "-E", "1e9", "--cpu", THREADS, "--tblout", tbl, "-o", os.devnull, hmm, allq],
        os.path.join(work, "hmmsearch.log"))
    best = {}
    for line in open(tbl):
        if line.startswith("#"):
            continue
        f = line.split()
        best[f[0]] = max(best.get(f[0], -1e9), float(f[5]))
    names = list(unaligned)
    sc = [best.get(n, 0.0) / lens[n] for n in names]
    zz = robust_z(sc)
    return {n for n, v in zip(names, zz) if v < -z}, dict(zip(names, sc)), bb
