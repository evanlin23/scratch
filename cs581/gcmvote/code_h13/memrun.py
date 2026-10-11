# AI-assisted (Claude), exploration code for CS581 project
"""Run a command; report wall, CPU, peak summed RSS of its process tree (sampled every 1 s) and the
largest single-process RSS (ru_maxrss of children).

    from memrun import run; run(cmd, log, cwd=None, env=None) -> dict
"""
import os
import resource
import subprocess
import threading
import time


def _tree_rss(root):
    kids = {}
    for p in os.listdir("/proc"):
        if not p.isdigit():
            continue
        try:
            with open("/proc/{}/stat".format(p)) as f:
                st = f.read().rsplit(")", 1)[1].split()
            kids.setdefault(int(st[1]), []).append(int(p))
        except OSError:
            continue
    total, todo = 0, [root]
    while todo:
        p = todo.pop()
        todo += kids.get(p, [])
        try:
            with open("/proc/{}/status".format(p)) as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        total += int(line.split()[1])
        except OSError:
            pass
    return total * 1024


def run(cmd, log, cwd=None, env=None):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.time()
    peak = [0]
    with open(log, "w") as f:
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, env=env)
        stop = threading.Event()

        def sample():
            while not stop.is_set():
                peak[0] = max(peak[0], _tree_rss(proc.pid))
                stop.wait(1.0)
        t = threading.Thread(target=sample, daemon=True)
        t.start()
        rc = proc.wait()
        stop.set()
        t.join()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    out = {"wall": round(time.time() - start, 1),
           "cpu": round(after.ru_utime - before.ru_utime + after.ru_stime - before.ru_stime, 1),
           "peak_tree_rss_gb": round(peak[0] / 2**30, 3), "max_proc_rss_gb": round(after.ru_maxrss / 2**20, 3),
           "rc": rc}
    if rc != 0:
        raise RuntimeError("failed ({}): {} -> {}".format(rc, " ".join(cmd), out))
    return out
