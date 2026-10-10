"""Aligner command lines (all from unaligned FASTA, 4 threads) and a timed runner.

Every tool writes FASTA to OUT. Binaries come from the micromamba env "bio"
(MUSCLE 5.3, FAMSA, TWILIGHT) or apt (MAFFT 7.505, Clustal Omega).
"""

import os
import signal
import subprocess
import time

BIO = "/opt/mm/root/envs/bio/bin"
MAFFT = "/usr/bin/mafft"
CLUSTALO = "/usr/bin/clustalo"
MUSCLE = "/opt/mm/root/envs/muscle5/bin/muscle"  # bio env's muscle is PASTA's bundled v3.8
FAMSA = os.path.join(BIO, "famsa")
TWI_ITER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "twilight_iter.py")


def _mafft(*flags):
    return lambda i, o, t, w: ([MAFFT, *flags, "--thread", str(t), i], o)  # stdout -> o


TOOLS = {
    # name: f(in, out, threads, workdir) -> (argv, stdout_path or None)
    "mafft-auto": _mafft("--auto"),
    "mafft-linsi": _mafft("--localpair", "--maxiterate", "1000"),
    "mafft-ginsi": _mafft("--globalpair", "--maxiterate", "1000"),
    "mafft-einsi": _mafft("--genafpair", "--ep", "0", "--maxiterate", "1000"),
    # MAGUS's own subset/backbone command (magus/tools/external_tools.py runMafft), system MAFFT 7.505
    "mafft-linsi-magus": _mafft("--localpair", "--maxiterate", "1000", "--ep", "0.123", "--quiet", "--anysymbol"),
    "mafft-parttree": _mafft("--retree", "2", "--parttree"),
    "mafft-fftns2": _mafft("--retree", "2"),
    "clustalo": lambda i, o, t, w: ([CLUSTALO, "-i", i, "-o", o, "--threads", str(t), "--force", "--outfmt", "fa"], None),
    "muscle5": lambda i, o, t, w: ([MUSCLE, "-align", i, "-output", o, "-threads", str(t)], None),
    "muscle5-super5": lambda i, o, t, w: ([MUSCLE, "-super5", i, "-output", o, "-threads", str(t)], None),
    "famsa": lambda i, o, t, w: ([FAMSA, "-t", str(t), i, o], None),
    "famsa-medoid": lambda i, o, t, w: ([FAMSA, "-t", str(t), "-medoidtree", i, o], None),
    "twilight": lambda i, o, t, w: (["python3", TWI_ITER, i, o, str(t), os.path.join(w, "twilight_tmp"), "3"], None),
    # diagnostics: same aligner, true tree as guide tree (WORKDIR/true_tree.nwk, written by diag runs)
    "twilight-truetree": lambda i, o, t, w: ([os.path.join(BIO, "twilight"), "-i", i, "-t", os.path.join(w, "true_tree.nwk"),
                                              "-o", o, "-C", str(t), "--overwrite"], None),
    "famsa-truetree": lambda i, o, t, w: ([FAMSA, "-t", str(t), "-gt", "import", os.path.join(w, "true_tree.nwk"), i, o], None),
    "twilight-1": lambda i, o, t, w: (["python3", TWI_ITER, i, o, str(t), os.path.join(w, "twilight_tmp"), "1"], None),
}


def run(tool, unaligned, out, threads, workdir, timeout, log, cwd=None):
    """Run TOOL; returns dict(status, wall, cpu, maxrss_mb)."""
    argv, stdout = TOOLS[tool](unaligned, out, threads, workdir)
    workdir = cwd or workdir
    start = time.time()
    with open(log, "w") as lf, open(stdout or os.devnull, "w") as so:
        proc = subprocess.Popen(argv, stdout=so if stdout else lf, stderr=lf, cwd=workdir, start_new_session=True)
        status = None
        while True:  # wait4 gives this child's own rusage (CPU, peak RSS of the largest descendant)
            pid, st, ru = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                rc = os.waitstatus_to_exitcode(st)
                status = status or ("ok" if rc == 0 else "fail rc={}".format(rc))
                break
            if status is None and time.time() - start > timeout:
                os.killpg(proc.pid, signal.SIGKILL)
                status = "timeout"
            time.sleep(0.05)
        proc.returncode = rc
    wall = time.time() - start
    if status == "ok" and (not os.path.exists(out) or os.path.getsize(out) == 0):
        status = "no output"
    return {"status": status, "wall": round(wall, 2), "cpu": round(ru.ru_utime + ru.ru_stime, 1),
            "maxrss_mb": round(ru.ru_maxrss / 1024)}
