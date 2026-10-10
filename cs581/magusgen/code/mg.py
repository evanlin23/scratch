"""One MAGUS variant for DNA, RNA and proteins? Merge-only paired comparisons vs MAGUS.

    python3 mg.py rep  NAME REP_DIR TRUE               # replicate from a cached MAGUS run (worker branches)
    python3 mg.py run  REP_DIR VARIANT [...]           # variants -> REP_DIR/results.jsonl

Builds on cs581/gcmgen/code/gg.py (and through it bbevidence/bbe.py and protcons/pc.py), whose
variant grammar is used unchanged for the one-stage merges:

  linsi                 MAGUS (control; reproduces MAGUS's own output)
  linsi&fftns2-op3      consensus evidence: only pairs L-INS-i and FFT-NS-2 --op 3 both make
  softW:A&B             soft consensus: confirmed pairs weight 1, unconfirmed W (gg.py)
  linsi|cons0.7         self-consistency mask (columns re-aligned by < 70% of the other backbones)

New here, two-stage "self-soft" (cs581/code/gcmx/e2e_bench.py) on top of any one-stage variant BASE:

  ss:BASE    BASE's merged alignment restricted to 8 random sequences per subset -> 10 pseudo-backbones,
             HMM-extended to all sequences (gcmx.extend); each subset alignment split into 3 similarity
             groups (gcmx.split); all 75 groups merged at once with the pseudo-backbones as the only
             evidence. ss:linsi is the published self-soft recipe.
  ssu:BASE   as ss:BASE, plus BASE's own (200-sequence) backbones as extra evidence.

Timing per row: bb_wall / prep_wall / merge_wall as in bbe.py; for ss rows, ev_wall (pseudo-backbones +
HMM extension), split_wall and merge_wall of the second stage (the first stage is the BASE row).
"""

import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gcmgen", "code"))
import gg  # noqa: E402  (imports and patches bbe, pc)
from gg import bbe  # noqa: E402
from gcmx import fasta  # noqa: E402
from gcmx.bbtool_bench import acc_ref  # noqa: E402
from gcmx.pilot import self_backbones  # noqa: E402

CODE = bbe.CODE
T = str(bbe.THREADS)


def results(rep):
    res = os.path.join(rep, "results.jsonl")
    rows = {}
    if os.path.exists(res):
        for l in open(res):
            r = json.loads(l)
            rows[r["variant"]] = r
    return rows


def base_out(rep, base):
    return os.path.join(rep, "variants", bbe.safe(base.replace(":", "_c_")), "out.fasta")


def sh(cmd, log):
    start = time.time()
    with open(log, "w") as f:
        subprocess.run(cmd, cwd=CODE, stdout=f, stderr=subprocess.STDOUT, check=True)
    return round(time.time() - start, 1)


def unaligned(rep):
    p = os.path.join(rep, "unaligned.fasta")
    if not os.path.exists(p):
        fasta.write(fasta.ungap(fasta.read(os.path.join(rep, "true.fasta"))), p + ".tmp")
        os.replace(p + ".tmp", p)
    return p


def split3(rep):
    d = os.path.join(rep, "split_m3")
    if not os.path.exists(d + ".DONE"):
        shutil.rmtree(d, ignore_errors=True)
        w = sh([sys.executable, "-m", "gcmx.split", os.path.join(rep, "inputs", "subalignments"), d, "3"],
               os.path.join(rep, "split_m3.log"))
        json.dump({"wall": w}, open(d + ".DONE", "w"))
    return d, json.load(open(d + ".DONE"))["wall"]


def self_soft(rep, name):
    kind, base = name.split(":", 1)
    rows = results(rep)
    if base not in rows or not os.path.exists(base_out(rep, base)):
        if base in rows:  # row exists but its output was cleaned: rerun the base
            drop(rep, base)
        gg.run(rep, [base])
    split, split_wall = split3(rep)
    vd = os.path.join(rep, "variants", bbe.safe(name.replace(":", "_c_")))
    shutil.rmtree(vd, ignore_errors=True)
    os.makedirs(vd)
    bb, ev = os.path.join(vd, "self_bb"), os.path.join(vd, "ev")
    start = time.time()
    self_backbones(base_out(rep, base), os.path.join(rep, "inputs", "subalignments"), bb)
    prep = time.time() - start
    ev_wall = round(prep + sh([sys.executable, "-m", "gcmx.extend", bb, unaligned(rep), ev, "--jobs", T],
                              os.path.join(vd, "extend.log")), 1)
    if kind == "ssu":
        files, _, _ = gg.files_of(rep, base)
        for lab, a in files:
            fasta.write(a, os.path.join(ev, "base_" + lab + ".txt"))
    out = os.path.join(vd, "out.fasta")
    m_wall = sh([sys.executable, "-m", "gcmx.run_magus", "-np", T, "--gcmx-mclthreads", T, "-d",
                 os.path.join(vd, "work"), "-s", split, "-b", ev, "-o", out], os.path.join(vd, "magus.log"))
    shutil.rmtree(os.path.join(vd, "work"), ignore_errors=True)
    shutil.rmtree(ev, ignore_errors=True)
    s = acc_ref(os.path.join(rep, "true.fasta"), out)
    b = results(rep)[base]
    row = {"rep": os.path.basename(rep.rstrip("/")), "variant": name, "base_merge_wall": b.get("merge_wall"),
           "base_prep_wall": b.get("prep_wall"), "base_bb_wall": b.get("bb_wall"), "base_bb_sum": b.get("bb_sum"),
           "ev_wall": ev_wall, "split_wall": split_wall, "merge_wall": m_wall, **s}
    with open(os.path.join(rep, "results.jsonl"), "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)


def drop(rep, name):
    res = os.path.join(rep, "results.jsonl")
    keep = [l for l in open(res) if json.loads(l)["variant"] != name]
    open(res + ".tmp", "w").writelines(keep)
    os.replace(res + ".tmp", res)


def run(rep, names):
    for name in names:
        if name in results(rep):
            continue
        if name.startswith(("ss:", "ssu:")):
            self_soft(rep, name)
        else:
            gg.run(rep, [name])


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "rep":
        gg.make_rep(sys.argv[2], os.path.abspath(sys.argv[3]), sys.argv[4])
    elif cmd == "run":
        run(os.path.abspath(sys.argv[2]), sys.argv[3:])
