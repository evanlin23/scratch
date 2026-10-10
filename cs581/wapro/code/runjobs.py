"""Run a job list (one JSON per line) with P parallel workers. Restartable (bench.py skips finished work).
Usage: python runjobs.py JOBS.txt P"""
import json, os, subprocess, sys
from multiprocessing.pool import ThreadPool
H = os.path.dirname(os.path.abspath(__file__))
PY = "/opt/mm/root/envs/gdl/bin/python"


def run(line):
    j = json.loads(line)
    if not os.path.exists(j["genes"]):
        return "missing " + j["genes"]
    subprocess.run([PY, os.path.join(H, "bench.py"), "--genes", j["genes"], "--true", j["true"], "--key", json.dumps(j["key"]),
                    "--out", j["out"], "--methods", j["methods"], "--threads", str(j["threads"])],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return "done " + json.dumps(j["key"])


if __name__ == "__main__":
    lines = [l for l in open(sys.argv[1]) if l.strip()]
    with ThreadPool(int(sys.argv[2])) as p:
        for i, r in enumerate(p.imap_unordered(run, lines)):
            print(i, r, flush=True)
