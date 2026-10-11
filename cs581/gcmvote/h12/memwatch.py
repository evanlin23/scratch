# AI-assisted (Claude), exploration code for CS581 project
"""Run a command, sample the summed RSS of its whole process tree every 2 s, enforce a wall-clock limit.

    python3 memwatch.py OUT.json LIMIT_SECONDS CMD [ARGS...]

Writes {"cmd", "wall", "peak_rss_gb", "rc", "timed_out"} to OUT.json (peak is of the sampled tree sum,
so it can miss spikes shorter than the 2 s interval)."""
import json
import os
import signal
import subprocess
import sys
import time


def tree_rss(root):
    kids = {}
    for p in os.listdir("/proc"):
        if not p.isdigit():
            continue
        try:
            with open("/proc/%s/stat" % p) as f:
                ppid = int(f.read().rsplit(")", 1)[1].split()[1])
            kids.setdefault(ppid, []).append(int(p))
        except (OSError, IndexError, ValueError):
            pass
    total, stack = 0, [root]
    while stack:
        p = stack.pop()
        stack += kids.get(p, [])
        try:
            with open("/proc/%d/status" % p) as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        total += int(line.split()[1]) * 1024
        except OSError:
            pass
    return total


def main():
    out, limit, cmd = sys.argv[1], float(sys.argv[2]), sys.argv[3:]
    start = time.time()
    proc = subprocess.Popen(cmd, start_new_session=True)
    peak, timed_out = 0, False
    while proc.poll() is None:
        peak = max(peak, tree_rss(proc.pid))
        if time.time() - start > limit:
            timed_out = True
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            break
        time.sleep(2)
    row = {"cmd": " ".join(cmd), "wall": round(time.time() - start, 1), "peak_rss_gb": round(peak / 2**30, 2),
           "rc": proc.returncode, "timed_out": timed_out}
    with open(out, "w") as f:
        json.dump(row, f)
    print(json.dumps(row), flush=True)
    sys.exit(0 if proc.returncode == 0 else 1)


if __name__ == "__main__":
    main()
