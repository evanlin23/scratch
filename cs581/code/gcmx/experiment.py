"""Run GCM merge variants on fixed inputs and record accuracy + time.

    python -m gcmx.experiment --dataset NAME --true TRUE.fa \
        --subalignments DIR --backbones DIR --outdir OUT \
        --variant count:"--gcmx-weight count" --variant frac:"--gcmx-weight frac" ...

Each variant reuses the same subalignments and backbones, so only the merge
step differs between variants. Results are appended to OUT/results.jsonl.
"""

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

from . import score


def run_variant(args, name, extra):
    workdir = os.path.join(args.outdir, args.dataset, name)
    output = os.path.join(args.outdir, args.dataset, name + ".fasta")
    # always start clean: MAGUS reuses graph/cluster files and waits on stale task
    # entries left in its working directory by an interrupted run
    shutil.rmtree(workdir, ignore_errors=True)
    if os.path.exists(output):
        os.remove(output)
    os.makedirs(workdir)
    extra_args = shlex.split(extra)
    cmd = [sys.executable, "-m", "gcmx.run_magus", "-np", str(args.threads), "-d", workdir, "-o", output]
    if "-s" not in extra_args:
        cmd += ["-s", args.subalignments]
    if args.backbones and "-b" not in extra_args:
        cmd += ["-b", args.backbones]
    cmd += extra_args
    start = time.time()
    with open(os.path.join(workdir, "stdout.log"), "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                              cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    elapsed = time.time() - start
    record = {"dataset": args.dataset, "variant": name, "args": extra, "seconds": round(elapsed, 1),
              "returncode": proc.returncode}
    log_path = os.path.join(workdir, "log.txt")
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                if "total cost of" in line:
                    record["cutCost"] = float(line.rsplit("total cost of", 1)[1])
                elif "aborted with an exception" in line:
                    record["aborted"] = True
    if os.path.exists(output):
        record.update(score.fastsp(args.true, output))
    # keep only logs: graph files of soft-constraint merges are ~1 GB each and fill the disk
    if not os.environ.get("GCMX_KEEP_WORKDIR"):
        for item in os.listdir(workdir):
            if item not in ("log.txt", "stdout.log"):
                path = os.path.join(workdir, item)
                shutil.rmtree(path) if os.path.isdir(path) else os.remove(path)
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--true", required=True)
    parser.add_argument("--subalignments", required=True, help="default -s; a variant may pass its own -s/-b")
    parser.add_argument("--backbones", default=None)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--variant", action="append", required=True, help="name:extra-args")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    variants = [v.split(":", 1) if ":" in v else (v, "") for v in args.variant]
    os.makedirs(args.outdir, exist_ok=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_variant, args, name, extra) for name, extra in variants]
        for future in futures:
            record = future.result()
            print(json.dumps(record), flush=True)
            with open(os.path.join(args.outdir, "results.jsonl"), "a") as f:
                f.write(json.dumps(record) + "\n")


if __name__ == "__main__":
    main()
