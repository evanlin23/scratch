"""Estimated alignment for Q3: WITCH (witch-msa 1.0.10) with a MAFFT --auto backbone (MAGUS optional), 1 CPU.

    python align_witch.py DATASET REP

Backbone = sequences within +-25% of the 3rd quartile of ungapped lengths (UPP's -M 0.75 rule; WITCH's
own default uses the median of ALL lengths, which on HF data falls between fragments and full-length
sequences), at most 1000 (random, seed 0). Backbone aligned with MAFFT --auto (WITCH_BACKBONE=magus: MAGUS); the other
sequences are added with WITCH (default: 10 HMMs, weighted). Insertion columns (lower-case in WITCH's
output) are removed, as in UPP's masked alignment.
Writes $MLDATA/DS/R<rep>/witch.fasta and witch.cost.json (cpu, wall, peak RSS of all steps).
"""
import json
import os
import random
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402

E = "/opt/mm/root/envs/aln"
PY = E + "/bin/python"
BACKBONE = os.environ.get("WITCH_BACKBONE", "mafft")  # amendment 3: MAGUS too slow on this VM
MAGUS = E + "/lib/python3.11/site-packages/witch_msa/tools/magus/magus.py"


def run(cmd, log, cost):
    t = time.time()
    with open(log, "w") as lf:
        p = subprocess.Popen(cmd, stdout=lf, stderr=lf, env=dict(os.environ, PATH=E + "/bin:" + os.environ["PATH"], OMP_NUM_THREADS="1"))
        _, st, ru = os.wait4(p.pid, 0)
    # RUSAGE of the direct child only covers grandchildren it waited for; take RUSAGE_CHILDREN deltas instead
    if os.waitstatus_to_exitcode(st) != 0:
        raise RuntimeError("failed: %s (%s)" % (cmd, log))
    cost["wall"] += time.time() - t
    cost["rss"] = max(cost["rss"], ru.ru_maxrss / 1024.0)


def main(ds, rep):
    d = os.path.join(rt.MLDATA, ds, "R%s" % rep)
    if os.path.exists(os.path.join(d, "witch.fasta")):
        return
    w = os.path.join(d, "witch_work")
    shutil.rmtree(w, ignore_errors=True)
    os.makedirs(w)
    names, seqs = rt.read_fasta(os.path.join(d, "true_align.fasta"))
    seqs = [s.upper().replace("-", "").replace(".", "") for s in seqs]
    L = sorted(len(s) for s in seqs)
    q3 = L[int(0.75 * len(L))]
    full = [i for i, s in enumerate(seqs) if 0.75 * q3 <= len(s) <= 1.25 * q3]
    if len(full) > 1000:
        full = sorted(random.Random(0).sample(full, 1000))
    fs = set(full)
    rt.write_fasta(os.path.join(w, "bb.fa"), [names[i] for i in full], [seqs[i] for i in full])
    rt.write_fasta(os.path.join(w, "q.fa"), [names[i] for i in range(len(names)) if i not in fs],
                   [seqs[i] for i in range(len(names)) if i not in fs])
    cost = {"cpu": 0.0, "wall": 0.0, "rss": 0.0}
    r0 = os.times()
    if BACKBONE == "magus":
        run([PY, MAGUS, "-i", os.path.join(w, "bb.fa"), "-o", os.path.join(w, "bb.aln"), "-d",
             os.path.join(w, "magus"), "-np", "1"], os.path.join(w, "magus.log"), cost)
    else:  # MAFFT --auto (FFT-NS-i on ~500 sequences), 1 thread
        t = time.time()
        with open(os.path.join(w, "bb.aln"), "w") as o, open(os.path.join(w, "mafft.log"), "w") as lf:
            p = subprocess.Popen(["/usr/bin/mafft", "--auto", "--thread", "1", "--quiet", os.path.join(w, "bb.fa")],
                                 stdout=o, stderr=lf)
            _, st, ru = os.wait4(p.pid, 0)
        assert os.waitstatus_to_exitcode(st) == 0
        cost["wall"] += time.time() - t
    run([PY, E + "/bin/witch.py", "-b", os.path.join(w, "bb.aln"), "-q", os.path.join(w, "q.fa"), "-d",
         os.path.join(w, "witch"), "-o", "aln.fasta", "-t", "1"], os.path.join(w, "witch.log"), cost)
    r1 = os.times()
    cost["cpu"] = (r1.children_user - r0.children_user) + (r1.children_system - r0.children_system)
    n2, s2 = rt.read_fasta(os.path.join(w, "witch", "aln.fasta"))
    keep = [j for j in range(len(s2[0])) if not any(s[j].islower() for s in s2)]
    rt.write_fasta(os.path.join(d, "witch.fasta.tmp"), n2, ["".join(s[j] for j in keep) for s in s2])
    os.replace(os.path.join(d, "witch.fasta.tmp"), os.path.join(d, "witch.fasta"))
    cost.update(backbone=len(full), q3=q3, cols=len(keep), cols_raw=len(s2[0]))
    json.dump(cost, open(os.path.join(d, "witch.cost.json"), "w"))
    print(ds, rep, cost)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
