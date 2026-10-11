# AI-assisted (Claude), exploration code for CS581 project
"""Helper h13, after magus_h13.py: baselines (gg.py), every pre-declared vote variant (gcmvote/code/run.py),
FastTree trees; each step under memrun (wall, CPU, peak memory). Restartable; rows go to OUT_DIR.

    python3 pipeline.py WORK OUT_DIR STEP [...]     STEP in: gg vote trees:VARIANT

  WORK/rep        gg.py replicate (gg.py writes REP/results.jsonl)
  WORK/rep_vote   the same inputs for run.py (separate dir: run.py's results.jsonl has a different schema)
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
from memrun import run  # noqa: E402

GG = ["linsi", "linsi#es3", "linsi#es4", "linsi#es5"]
VOTE = ["magus", "es4", "hard", "soft", "soft2", "soft4", "hard-bb", "soft-bb",
        "frac0.2", "frac0.3", "frac0.4", "frac0.5", "hard+mask"]
CODE = os.path.join(REPO, "code")
ENV = dict(os.environ, PYTHONPATH=CODE)


def done_of(path, key):
    return {json.loads(l)[key] for l in open(path)} if os.path.exists(path) else set()


def append(path, row):
    with open(path, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)


def main(work, outdir, steps):
    rep, rv = os.path.join(work, "rep"), os.path.join(work, "rep_vote")
    os.makedirs(outdir, exist_ok=True)
    aln_out = os.path.join(outdir, "aln.jsonl")
    magus = json.load(open(os.path.join(work, "d0", "state.json")))
    base = {"dataset": "RNASim10k_R0", "nseq": magus["nseq"], "helper": "h13"}
    done = {(r["source"], r["variant"]) for r in map(json.loads, open(aln_out))} if os.path.exists(aln_out) else set()
    if ("magus_e2e", "magus") not in done:
        m = magus["magus"]
        append(aln_out, {**base, "source": "magus_e2e", "variant": "magus", **m,
                         "avgErr_pct": round(100 * (m["SPFN"] + m["SPFP"]) / 2, 3)})
    if "merge-mafft-slow" in magus and ("magus_merge_slowgraph", "magus") not in done:
        m = magus["merge-mafft-slow"]
        append(aln_out, {**base, "source": "magus_merge_slowgraph", "variant": "magus", **m,
                         "avgErr_pct": round(100 * (m["SPFN"] + m["SPFP"]) / 2, 3)})
    for step in steps:
        if step == "gg":
            for v in GG:
                if ("gg", v) in done:
                    continue
                log = os.path.join(work, "gg_{}.log".format(v.replace("#", "_")))
                try:
                    m = run([sys.executable, os.path.join(REPO, "gcmgen", "code", "gg.py"), "run", rep, v], log,
                            cwd=CODE, env=ENV)
                except RuntimeError as e:
                    append(aln_out, {**base, "source": "gg", "variant": v, "failed": str(e)[-500:]})
                    continue
                r = [x for x in map(json.loads, open(os.path.join(rep, "results.jsonl"))) if x["variant"] == v][-1]
                append(aln_out, {**base, "source": "gg", **r, **m,
                                 "avgErr_pct": round(100 * (r["SPFN"] + r["SPFP"]) / 2, 3)})
        elif step == "vote":
            for sub in ("inputs", "true.fasta"):
                os.makedirs(rv, exist_ok=True)
                if not os.path.exists(os.path.join(rv, sub)):
                    os.symlink(os.path.join(rep, sub), os.path.join(rv, sub))
            for v in VOTE:
                if ("vote", v) in done:
                    continue
                log = os.path.join(work, "vote_{}.log".format(v.replace("+", "_")))
                try:
                    m = run([sys.executable, os.path.join(REPO, "gcmvote", "code", "run.py"), rv, v], log,
                            cwd=CODE, env=ENV)
                except RuntimeError as e:
                    append(aln_out, {**base, "source": "vote", "variant": v, "failed": str(e)[-500:]})
                    continue
                rows = [x for x in map(json.loads, open(os.path.join(rv, "results.jsonl"))) if x["variant"] == v] \
                    if os.path.exists(os.path.join(rv, "results.jsonl")) else []
                if not rows:
                    append(aln_out, {**base, "source": "vote", "variant": v, **m,
                                     "failed": open(log).read()[-1500:]})
                    continue
                r = rows[-1]
                r.pop("cutoff_by_n", None)
                append(aln_out, {**base, "source": "vote", **r, **m,
                                 "cutoff_by_n": json.load(open(os.path.join(rv, "vote", v, "model.json")))["cutoff_by_n"],
                                 "avgErr_pct": round(100 * (r["SPFN"] + r["SPFP"]) / 2, 3)})
        elif step.startswith("trees"):
            alns = {"true": os.path.join(rep, "true.fasta"), "magus": os.path.join(work, "d0", "magus.fasta"),
                    "es4": os.path.join(rep, "variants", "linsi_es_4", "out.fasta")}
            for v in step.split(":")[1:]:
                alns["vote-" + v] = os.path.join(rv, "vote", v, "out.fasta")
            for meth, path in alns.items():
                if not os.path.exists(path):
                    print("missing", meth, path, flush=True)
                    continue
                tlog = os.path.join(work, "tree_{}.log".format(meth))
                m = run(["/opt/mm/root/envs/pasta183/bin/python", os.path.join(HERE, "trees.py"), "RNASim10k_R0",
                         "/opt/data/Datasets/RNASim/10000/R0/true_tree.tre", os.path.join(outdir, "trees.jsonl"),
                         "{}={}".format(meth, path)], tlog)
                print(json.dumps({"tree": meth, **m}), flush=True)


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3:])
