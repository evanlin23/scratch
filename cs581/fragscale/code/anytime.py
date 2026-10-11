"""Run a tree search and record its anytime trajectory (best tree so far vs CPU time).

    snaps = run_anytime(cmd, log, snapshot, snapdir, cap_cpu=None, cpu_offset=0.0)

`snapshot()` returns the current best tree as a Newick string (or None); it is polled every
POLL seconds and every new tree is saved as snapdir/NNNN.tre with the process' CPU time
(user+sys from /proc, plus cpu_offset for resumed runs) in snapdir/index.tsv:
    cpu_s  wall_s  file
If cap_cpu is given, the process is killed (SIGTERM) once its CPU exceeds the cap.
Returns (Cost(cpu, wall, rss), killed).

Snapshot sources:
  rx_snapshot(prefix): RAxML-NG rewrites <prefix>.raxml.lastTree.TMP at every checkpoint
  iq_snapshot(prefix): IQ-TREE checkpoint <prefix>.ckp.gz, best tree of the candidate set
"""
import gzip
import os
import re
import signal
import subprocess
import time

POLL = 2.0
TCK = os.sysconf("SC_CLK_TCK")


class Cost:
    def __init__(self, cpu=0.0, wall=0.0, rss=0.0):
        self.cpu, self.wall, self.rss = cpu, wall, rss

    def add(self, o):
        return Cost(self.cpu + o.cpu, self.wall + o.wall, max(self.rss, o.rss))

    def d(self):
        return {"cpu": round(self.cpu, 1), "wall": round(self.wall, 1), "rss": round(self.rss, 1)}


def proc_cpu(pid):
    try:
        f = open("/proc/%d/stat" % pid).read()
    except OSError:
        return None
    x = f[f.rindex(")") + 2:].split()
    return (int(x[11]) + int(x[12])) / TCK  # utime + stime


def rx_snapshot(prefix):
    seen = {"last": False}

    def f():
        # the parsimony start tree only until the first checkpoint tree exists (RAxML-NG deletes
        # lastTree.TMP when it finishes; the caller then records the final tree itself)
        for suf in (".raxml.lastTree.TMP",) + (() if seen["last"] else (".raxml.startTree",)):
            p = prefix + suf
            if os.path.exists(p):
                try:
                    s = open(p).read().strip()
                except OSError:
                    continue
                if s.endswith(";"):
                    seen["last"] = seen["last"] or suf == ".raxml.lastTree.TMP"
                    return s
        return None
    return f


def iq_snapshot(prefix):
    """Best-scoring tree in IQ-TREE's checkpoint (CandidateSet entries are 'score tree')."""
    def f():
        p = prefix + ".ckp.gz"
        if not os.path.exists(p):
            return None
        try:
            txt = gzip.open(p, "rt").read()
        except (OSError, EOFError):
            return None
        best = None
        for m in re.finditer(r'^\s+\d+: "?(-[0-9.]+) (\(.*;)"?\s*$', txt, re.M):
            sc = float(m.group(1))
            if best is None or sc > best[0]:
                best = (sc, m.group(2))
        if best:
            return best[1]
        m = re.search(r'^initTree: (\(.*;)\s*$', txt, re.M)  # before the candidate set exists
        return m.group(1) if m else None
    return f


def run_anytime(cmd, log, snapshot, snapdir, cap_cpu=None, cpu_offset=0.0, cwd=None):
    os.makedirs(snapdir, exist_ok=True)
    idx = os.path.join(snapdir, "index.tsv")
    n = sum(1 for _ in open(idx)) if os.path.exists(idx) else 0
    last = None
    t0 = time.time()
    killed = False
    with open(log, "a") as lf:
        p = subprocess.Popen(cmd, stdout=lf, stderr=lf, cwd=cwd)
        cpu = 0.0
        while True:
            time.sleep(POLL)
            c = proc_cpu(p.pid)
            if c is not None:
                cpu = c
            s = snapshot()
            if s and s != last:
                last = s
                fn = os.path.join(snapdir, "%04d.tre" % n)
                open(fn, "w").write(s + "\n")
                with open(idx, "a") as f:
                    f.write("%.1f\t%.1f\t%s\n" % (cpu + cpu_offset, time.time() - t0, os.path.basename(fn)))
                n += 1
            pid, status, ru = os.wait4(p.pid, os.WNOHANG)
            if pid != 0:
                break
            if cap_cpu is not None and cpu + cpu_offset > cap_cpu and not killed:
                p.send_signal(signal.SIGTERM)
                killed = True
    p.returncode = os.waitstatus_to_exitcode(status)
    return Cost(ru.ru_utime + ru.ru_stime, time.time() - t0, ru.ru_maxrss / 1024.0), killed, p.returncode


def at(snapdir, T):
    """Path of the last snapshot with cpu <= T, or None."""
    best = None
    idx = os.path.join(snapdir, "index.tsv")
    if not os.path.exists(idx):
        return None
    for line in open(idx):
        c, w, f = line.split()
        if float(c) <= T:
            best = os.path.join(snapdir, f)
    return best
