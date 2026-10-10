"""PASTA end to end: default merger (OPAL) vs a pairmerge merger plugged in via the MUSCLE slot.

    python pasta_e2e.py RESULTS.jsonl ARM[,ARM] DATASET:REP_DIR ...   (ARM: opal | progdp | gcm)

Both arms use the same PASTA config (PASTA 1.9 defaults for DNA: centroid decomposition,
max subproblem 200, MAFFT L-INS-i subsets, FastTree) except the merger; `--iter-limit`
from PASTA_ITERS (default 1, budget). The progdp/gcm arms run cs581/pairmerge/code/pasta_shim.sh
as PASTA's "muscle" with backbones of PAIRMERGE_R x PAIRMERGE_M sequences per merge.
"""

import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))
from gcmx import fasta, score  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.environ.get("PAIRMERGE_RUNS", "/opt/runs/pairmerge_pasta")
BIO = "/opt/mm/root/envs/bio/bin"
TOOLS = BIO + "/sate-tools-linux/"


def config(path, merger_shim):
    """PASTA's exported default config with tool paths fixed (bioconda layout) and optionally the shim."""
    src = os.path.join(HERE, "pasta_default_config.txt")
    text = open(src).read().replace("/opt/mm/root/envs/bio/lib/python3.7/site-packages/bin/", TOOLS)
    if merger_shim:
        text = text.replace("path = " + TOOLS + "muscle\n", "path = " + os.path.join(HERE, "pasta_shim.sh") + "\n")
    open(path, "w").write(text)


def main():
    results, arms, reps = sys.argv[1], sys.argv[2].split(","), sys.argv[3:]
    done = set()
    if os.path.exists(results):
        done = {(r["dataset"], r["arm"]) for r in map(json.loads, open(results)) if "avgErr" in r}
    iters = os.environ.get("PASTA_ITERS", "1")
    for rep in reps:
        dataset, rep_dir = rep.split(":", 1)
        aln = next(os.path.join(rep_dir, f) for f in ("rose.aln.true.fasta", "true_align.txt")
                   if os.path.exists(os.path.join(rep_dir, f)))
        for arm in arms:
            if (dataset, arm) in done:
                continue
            work = os.path.join(RUNS, dataset, arm)
            os.makedirs(work, exist_ok=True)
            true = fasta.upper(fasta.read(aln))
            sub = int(os.environ.get("PASTA_SUBSAMPLE", "0"))  # budget: random n-taxon subsample (same for both arms)
            if sub:
                import random
                true = fasta.restrict(true, sorted(random.Random(dataset).sample(sorted(true), sub)))
            fasta.write(true, os.path.join(work, "true.fa"))
            fasta.write(fasta.ungap(true), os.path.join(work, "in.fa"))
            cfg = os.path.join(work, "cfg.txt")
            config(cfg, arm != "opal")
            env = dict(os.environ, PATH=BIO + ":" + os.environ["PATH"], PAIRMERGE_METHOD=arm,
                       PAIRMERGE_TMP=work, PAIRMERGE_LOG=os.path.join(work, "merges.log"))
            out = os.path.join(work, "out_{}".format(int(time.time())))
            cmd = [BIO + "/run_pasta.py", cfg, "-i", os.path.join(work, "in.fa"), "-d", "dna", "-o", out, "-j", "p",
                   "--iter-limit", iters, "--num-cpus", os.environ.get("PASTA_CPUS", "4"),
                   "--merger", "opal" if arm == "opal" else "muscle"]
            start = time.time()
            with open(os.path.join(work, "pasta.log"), "w") as log:
                proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=work)
            rec = {"dataset": dataset, "arm": arm, "iters": int(iters), "ntaxa": len(true), "seconds": round(time.time() - start, 1),
                   "returncode": proc.returncode}
            final = os.path.join(out, "p.marker001.in.fa.aln")
            if os.path.exists(final):
                rec.update(score.fastsp(os.path.join(work, "true.fa"), final))
            print(json.dumps(rec), flush=True)
            with open(results, "a") as f:
                f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
