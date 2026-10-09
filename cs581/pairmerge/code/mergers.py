"""Two-alignment mergers behind one interface.

    merge(name, a_fasta, b_fasta, out_fasta, workdir, unaligned=None, threads=1, backbones=None)

Every merger takes two fixed alignments A and B (disjoint taxa) and writes one
alignment of A u B. The columns of A and B are kept intact (checked by
`check_constraints`): only the A-B column matching and the gap layout differ.

  mafft-merge   MAFFT --merge (sub-MSAs kept fixed, default FFT-NS-2 guide)
  mafft-merge-l MAFFT --localpair --merge (G-INS-i style pairwise scores across sub-MSAs)
  muscle3       MUSCLE 3.8 -profile (as PASTA calls it)
  opal          OPAL profile-profile (as PASTA calls it; PASTA's default merger)
  gcm           MAGUS GCM on the two subsets (MCL + minclusters trace), backbones built by MAGUS
  progdp        exact two-alignment MWT DP (cs581/code/gcmx/progressive.py) on the same graph
"""

import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GCMX_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "code"))  # contains the gcmx package
sys.path.insert(0, GCMX_ROOT)
from gcmx import fasta  # noqa: E402

BIO = "/opt/mm/root/envs/bio"
PASTA_BIN = os.environ.get("PASTA_TOOLS", "")  # set by find_pasta_tools()
MAGUS_FLAGS = ("--graphbuildmethod mafft --graphbuildhmmextend false --graphclustermethod mcl "
               "--graphtracemethod minclusters --graphtraceoptimize false -r 10 -m 200 -f 4").split()


def find_pasta_tools():
    """PASTA ships muscle (3.8) and opal.jar in its bundled tools directory."""
    hits = glob.glob(os.path.join(BIO, "**", "opal.jar"), recursive=True)
    opal = hits[0] if hits else None
    muscle = None
    for cand in glob.glob(os.path.join(BIO, "**", "muscle*"), recursive=True):
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            v = subprocess.run([cand, "-version"], capture_output=True, text=True)
            if "v3.8" in (v.stdout + v.stderr):
                muscle = cand
                break
    return opal, muscle


OPAL, MUSCLE3 = find_pasta_tools()


def _run(cmd, log, **kw):
    with open(log, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        subprocess.run(cmd, stdout=kw.pop("stdout", f), stderr=f, check=True, **kw)


def mafft_merge(a, b, out, work, threads, local=False):
    A, B = fasta.read(a), fasta.read(b)
    both = os.path.join(work, "both.fa")
    fasta.write({**A, **B}, both)
    table = os.path.join(work, "table")
    with open(table, "w") as f:
        f.write(" ".join(str(i + 1) for i in range(len(A))) + "\n")
        f.write(" ".join(str(len(A) + i + 1) for i in range(len(B))) + "\n")
    cmd = ["mafft", "--thread", str(threads), "--quiet"]
    if local:
        cmd += ["--localpair", "--maxiterate", "0"]
    cmd += ["--merge", table, both]
    with open(out, "w") as o:
        _run(cmd, os.path.join(work, "log.txt"), stdout=o)


def muscle3(a, b, out, work, threads):
    _run([MUSCLE3, "-profile", "-in1", a, "-in2", b, "-out", out, "-quiet"], os.path.join(work, "log.txt"))


def opal(a, b, out, work, threads):
    # PASTA's call: java -Xmx.. -jar opal.jar --in A --in2 B --out OUT --align_method profile
    _run(["java", "-Xmx4g", "-jar", OPAL, "--in", a, "--in2", b, "--out", out, "--align_method", "profile"],
         os.path.join(work, "log.txt"))


def _subset_dir(a, b, work):
    d = os.path.join(work, "subsets")
    os.makedirs(d, exist_ok=True)
    shutil.copy(a, os.path.join(d, "subalignment_subset_1.txt"))
    shutil.copy(b, os.path.join(d, "subalignment_subset_2.txt"))
    return d


def gcm(a, b, out, work, threads, backbones=None, trace="minclusters", keep_backbones=None):
    """MAGUS GCM on the two subsets. backbones=None -> MAGUS builds them (-r 10 -m 200, MAFFT)."""
    subs = _subset_dir(a, b, work)
    mwork = os.path.join(work, "magus")
    shutil.rmtree(mwork, ignore_errors=True)
    flags = list(MAGUS_FLAGS)
    if trace == "progdp":
        flags[flags.index("mcl")] = "none"
        flags[flags.index("minclusters")] = "progdp"
    cmd = [sys.executable, "-m", "gcmx.run_magus", "-np", str(threads), "-d", mwork, "-s", subs, "-o", out] + flags
    if backbones:
        cmd += ["-b", backbones]
    _run(cmd, os.path.join(work, "log.txt"), cwd=GCMX_ROOT)
    if keep_backbones:
        os.makedirs(keep_backbones, exist_ok=True)
        for p in glob.glob(os.path.join(mwork, "graph", "backbone_*_mafft.txt")):
            shutil.copy(p, keep_backbones)
    shutil.rmtree(mwork, ignore_errors=True)


def check_constraints(a, b, out):
    """True iff the merged alignment induces exactly A and B (columns kept intact)."""
    M = fasta.upper(fasta.read(out))
    for path in (a, b):
        S = fasta.upper(fasta.read(path))
        if set(S) - set(M):
            return False
        induced = fasta.restrict(M, list(S))
        if any(induced[t] != S[t] for t in S):
            return False
    return len(M) == len(fasta.read(a)) + len(fasta.read(b))
