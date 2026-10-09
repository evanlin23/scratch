"""Pilot of learned edge weights (gcmx.learnweights) on the fan-out replicates.

    python -m gcmx.lw_pilot WORKDIR [--train R0] [--threads 2]

For every replicate finished by a fan-out worker (cached MAGUS subset
alignments + 10 MAFFT backbones in cs581/experiments/runs/NAME/inputs.tar.xz on
the claude/cs581-worker-* branches):
  1. default merge (MAGUS's count weights) while dumping edge features + labels
  2. train the model on the training replicates only (default: every *_R0 of the
     ROSE and RNASim conditions; BAliBASE and 16S are never trained on)
  3. merge every held-out replicate with learned weights
Results (SPFN/SPFP, seconds) go to WORKDIR/results.jsonl; one row per run.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time

from . import fasta, score

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(os.path.dirname(CODE))


def git(*args):
    return subprocess.run(["git", "-C", REPO] + list(args), capture_output=True, check=False)


def replicates():
    """name -> (worker ref, true alignment path)"""
    truth = {}
    for path in glob.glob(os.path.join(CODE, "fanout", "jobs_w*.txt")):
        for line in open(path):
            if line.strip():
                name, src = line.split()[:2]
                truth[name] = src if os.path.isabs(src) else os.path.join(REPO, src)
    git("fetch", "-q", "origin", "+refs/heads/claude/cs581-worker-*:refs/remotes/origin/claude/cs581-worker-*")
    refs = git("for-each-ref", "--format=%(refname:short)", "refs/remotes/origin/claude/cs581-worker-*").stdout.decode().split()
    found = {}
    for ref in refs:
        out = git("ls-tree", "--name-only", ref + ":cs581/experiments/runs").stdout.decode().split()
        for name in out:
            if name in truth and git("cat-file", "-e", "{}:cs581/experiments/runs/{}/inputs.tar.xz".format(ref, name)).returncode == 0:
                found[name] = (ref, truth[name])
    return found


def prepare(work, name, ref, true_src):
    rep = os.path.join(work, name)
    if os.path.exists(os.path.join(rep, "true.fasta")):
        return rep
    os.makedirs(rep, exist_ok=True)
    blob = git("show", "{}:cs581/experiments/runs/{}/inputs.tar.xz".format(ref, name)).stdout
    tar = os.path.join(rep, "inputs.tar.xz")
    open(tar, "wb").write(blob)
    subprocess.run(["tar", "xJf", tar, "-C", rep], check=True)
    fasta.write(fasta.upper(fasta.read(true_src)), os.path.join(rep, "true.fasta"))
    return rep


def merge(rep, out, threads, env_extra):
    work = out[:-6] + "_work"
    shutil.rmtree(work, ignore_errors=True)
    env = dict(os.environ, **env_extra)
    start = time.time()
    with open(out[:-6] + ".log", "w") as log:
        subprocess.run([sys.executable, "-m", "gcmx.run_magus", "-np", str(threads), "-d", work,
                        "-s", os.path.join(rep, "inputs", "subalignments"), "-b", os.path.join(rep, "inputs", "backbones"),
                        "-o", out], cwd=CODE, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    seconds = round(time.time() - start, 1)
    shutil.rmtree(work, ignore_errors=True)
    s = score.fastsp(os.path.join(rep, "true.fasta"), out)
    return {"seconds": seconds, **{k: s[k] for k in ("SPFN", "SPFP", "avgErr", "LenEst", "LenRef")}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("work")
    parser.add_argument("--train", default="R0")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--model-name", default="model")
    args = parser.parse_args()
    os.makedirs(args.work, exist_ok=True)
    results = os.path.join(args.work, "results.jsonl")
    done = set()
    if os.path.exists(results):
        done = {(r["dataset"], r["variant"]) for r in map(json.loads, open(results))}

    def record(name, variant, row):
        with open(results, "a") as f:
            f.write(json.dumps({"dataset": name, "variant": variant, **row}) + "\n")
        print(name, variant, json.dumps(row), flush=True)

    reps = replicates()
    is_train = lambda n: n.endswith("_" + args.train) and not n.startswith(("BBA", "16S"))
    print("replicates:", len(reps), "train:", sorted(n for n in reps if is_train(n)), flush=True)

    # 1. default merges + training dumps
    for name, (ref, true_src) in sorted(reps.items()):
        rep = prepare(args.work, name, ref, true_src)
        dump = os.path.join(rep, "edges.npz")
        if (name, "default") not in done or not os.path.exists(dump):
            row = merge(rep, os.path.join(rep, "default.fasta"), args.threads,
                        {"GCMX_EDGE_DUMP": dump, "GCMX_EDGE_TRUE": os.path.join(rep, "true.fasta")})
            if (name, "default") not in done:
                record(name, "default", row)

    # 2. train
    model = os.path.join(args.work, args.model_name + ".joblib")
    train = sorted(os.path.join(args.work, n, "edges.npz") for n in reps if is_train(n))
    if not os.path.exists(model):
        subprocess.run([sys.executable, "-m", "gcmx.learnweights", "train", model] + train, cwd=CODE, check=True)

    # 3. held-out merges with learned weights
    variant = "learned-" + args.model_name
    for name in sorted(reps):
        if is_train(name) or (name, variant) in done:
            continue
        rep = os.path.join(args.work, name)
        record(name, variant, merge(rep, os.path.join(rep, variant + ".fasta"), args.threads, {"GCMX_EDGE_MODEL": model}))


if __name__ == "__main__":
    main()
