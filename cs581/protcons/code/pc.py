"""Held-out test of consistency-filtered GCM evidence (see ../PREREG.md).

    python3 pc.py rep  WORK/NAME_d0 REP_DIR       # bbe.py replicate from a gcmx.bbtool_bench MAGUS run
    python3 pc.py run  REP_DIR VARIANT [...]      # merge-only variants -> REP_DIR/results.jsonl
    python3 pc.py gate REP_DIR                    # reference-free support statistic -> REP_DIR/gate.json
    python3 pc.py check REP_DIR                   # gate == pairdiff.py's support_a_only (needs full reference)

Reuses cs581/bbevidence/code/bbe.py (variants, masking, intersection, merge) unchanged. Differences from
bbe.run: the score is computed on the estimate restricted to the reference's sequences (HomFam: the seeds),
and the cross-subset scores only when the reference covers every sequence.
"""

import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bbevidence", "code"))
import bbe  # noqa: E402
import pairdiff  # noqa: E402
from gcmx import fasta  # noqa: E402
from gcmx.bbtool_bench import acc_ref  # noqa: E402

TAU = 0.615


def make_rep(work, rep):
    """REP_DIR/inputs/{subalignments,backbones (L-INS-i, *_mafft.txt only)}, true.fasta, sets/s0."""
    if os.path.exists(os.path.join(rep, "sets", "s0")):
        return
    shutil.rmtree(rep, ignore_errors=True)
    os.makedirs(os.path.join(rep, "inputs", "backbones"))
    shutil.copytree(os.path.join(work, "inputs", "subalignments"), os.path.join(rep, "inputs", "subalignments"))
    for f in os.listdir(os.path.join(work, "inputs", "backbones")):
        if f.endswith("_mafft.txt"):
            shutil.copy(os.path.join(work, "inputs", "backbones", f), os.path.join(rep, "inputs", "backbones"))
    bbe.prep(rep, os.path.join(work, "true.fasta"))
    shutil.copy(os.path.join(work, "unaligned.fasta"), os.path.join(rep, "unaligned.fasta"))
    st = json.load(open(os.path.join(work, "state.json")))
    json.dump({"magus": st["magus"], "nseq": st["nseq"], "nref": st["nref"]}, open(os.path.join(rep, "magus.json"), "w"))


def full_ref(rep):
    ref = fasta.read(os.path.join(rep, "true.fasta"))
    n = sum(len(a) for _, a in bbe.subsets(rep))
    return len(ref) == n


def run(rep, names):
    res = os.path.join(rep, "results.jsonl")
    done = set()
    if os.path.exists(res):
        done = {json.loads(l)["variant"] for l in open(res)}
    full = full_ref(rep)
    for name in names:
        if name in done:
            continue
        start = time.time()
        base, k = (name.split("#es") + [None])[:2]
        files, bb_wall, bb_sum = bbe.parse_variant(rep, base)
        prep_wall = round(time.time() - start, 1)
        out, m_wall = merge_edgesup(rep, name, files, int(k)) if k else bbe.merge(rep, name, files)
        s = acc_ref(os.path.join(rep, "true.fasta"), out)
        row = {"rep": os.path.basename(rep.rstrip("/")), "variant": name, "nbb": len(files), "bb_wall": bb_wall,
               "bb_sum": bb_sum, "prep_wall": prep_wall, "merge_wall": m_wall, **s}
        if full:
            row.update(bbe.cross_scores(rep, out))
        with open(res, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


def merge_edgesup(rep, name, files, k):
    """bbe.merge, but the GCM graph keeps only edges supported by >= k backbones (run_edgesup.py)."""
    vd = os.path.join(rep, "variants", bbe.safe(name).replace("#", "_es_"))
    shutil.rmtree(vd, ignore_errors=True)
    bb = os.path.join(vd, "bb")
    os.makedirs(bb)
    for lab, a in files:
        fasta.write(a, os.path.join(bb, lab + ".txt"))
    out = os.path.join(vd, "out.fasta")
    start = time.time()
    with open(os.path.join(vd, "magus.log"), "w") as log:
        subprocess.run([sys.executable, os.path.join(HERE, "run_edgesup.py"), str(k), "--gcmx-fastgraph", "false",
                        "-np", str(bbe.THREADS), "-d", os.path.join(vd, "work"),
                        "-s", os.path.join(rep, "inputs", "subalignments"), "-b", bb, "-o", out] + bbe.MERGE_FLAGS,
                       cwd=bbe.CODE, stdout=log, stderr=subprocess.STDOUT, check=True)
    wall = round(time.time() - start, 1)
    shutil.rmtree(os.path.join(vd, "work"), ignore_errors=True)
    return out, wall


def new_sets(rep, seed, n=10, per_subset=8):
    """bbe.new_sets, but sequences come from rep/unaligned.fasta (HomFam: true.fasta holds only the seeds)."""
    import random
    d = os.path.join(rep, "sets", "s{}".format(seed))
    if os.path.isdir(d) and len(os.listdir(d)) == n:
        return d
    os.makedirs(d, exist_ok=True)
    src = os.path.join(rep, "unaligned.fasta")
    unal = fasta.ungap(fasta.upper(fasta.read(src if os.path.exists(src) else os.path.join(rep, "true.fasta"))))
    rng = random.Random(1000 + seed)
    subs = [list(s) for _, s in bbe.subsets(rep)]
    for b in range(n):
        taxa = [t for s in subs for t in rng.sample(s, min(per_subset, len(s)))]
        fasta.write({t: unal[t] for t in taxa}, os.path.join(d, "backbone_{}.fa".format(b + 1)))
    return d


bbe.new_sets = new_sets


class FreeRep:
    """What pairdiff.pairs_of needs, without a reference: subset of every taxon (dummy reference columns)."""

    def __init__(self, rep):
        self.subset_of, self.refcol = {}, {}
        for i, (_, a) in enumerate(bbe.subsets(rep)):
            for t, s in a.items():
                self.subset_of[t] = i
                self.refcol[t] = np.zeros(sum(c not in "-." for c in s), dtype=np.int64)
        self.ref_gapfrac = np.zeros(1)


def support(R, A, B):
    """pairdiff.main's support_a_only / support_shared, verbatim logic."""
    gid, off = {}, 0
    for t, rc in R.refcol.items():
        gid[t] = np.arange(off, off + len(rc), dtype=np.int64)
        off += len(rc)
    G = off
    colmaps = []
    for _, a in A:
        cm = np.full(G, -1, dtype=np.int64)
        rc = bbe.residue_columns(a)
        for t in a:
            cm[gid[t]] = rc[t]
        colmaps.append(cm)
    sup = {"a_only": [0, 0], "shared": [0, 0]}
    n_a = n_only = 0
    for bi, ((_, a), (_, b)) in enumerate(zip(A, B)):
        ka, _, _ = pairdiff.pairs_of(R, a, gid)
        kb, _, _ = pairdiff.pairs_of(R, b, gid)
        inb = np.isin(ka, kb, assume_unique=True)
        n_a += len(ka)
        n_only += int((~inb).sum())
        x, y = ka // G, ka % G
        for name, sel in (("a_only", ~inb), ("shared", inb)):
            xs, ys = x[sel], y[sel]
            for bj, cm in enumerate(colmaps):
                if bj == bi:
                    continue
                cx, cy = cm[xs], cm[ys]
                both = (cx >= 0) & (cy >= 0)
                sup[name][0] += int((cx[both] == cy[both]).sum())
                sup[name][1] += int(both.sum())
    return {"support_a_only": round(sup["a_only"][0] / max(sup["a_only"][1], 1), 4),
            "support_shared": round(sup["shared"][0] / max(sup["shared"][1], 1), 4),
            "frac_a_only": round(n_only / max(n_a, 1), 4)}


def gate(rep):
    path = os.path.join(rep, "gate.json")
    if os.path.exists(path):
        return json.load(open(path))
    start = time.time()
    A, _, _ = bbe.parse_variant(rep, "linsi")
    B, bw, bs = bbe.parse_variant(rep, "clustalo")
    g = support(FreeRep(rep), A, B)
    g.update({"tau": TAU, "filter": g["support_a_only"] < TAU, "gate_wall": round(time.time() - start, 1),
              "clustalo_bb_wall": bw, "clustalo_bb_sum": bs})
    json.dump(g, open(path, "w"))
    print(json.dumps(g), flush=True)
    return g


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "rep":
        make_rep(os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]))
    elif cmd == "run":
        run(os.path.abspath(sys.argv[2]), sys.argv[3:])
    elif cmd == "gate":
        gate(os.path.abspath(sys.argv[2]))
    elif cmd == "check":
        rep = os.path.abspath(sys.argv[2])
        print("reference-free:", gate(rep))
        pairdiff.main(rep, "linsi", "clustalo")
