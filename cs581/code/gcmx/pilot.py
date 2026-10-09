"""Pilot study: run all merge variants + oracle conditions on prepared replicates.

    python -m gcmx.pilot RESULTS_DIR REPLICATE_DIR [REPLICATE_DIR ...]

Each REPLICATE_DIR comes from gcmx.prep (true.fasta, inputs/subalignments,
inputs/backbones). (dataset, variant) pairs already in
RESULTS_DIR/results.jsonl are skipped, so the script can be rerun as more
replicates finish. PILOT_ONLY=name1,name2 restricts the variants (for large
datasets); PILOT_JOBS sets how many variants run concurrently.
"""

import glob
import json
import os
import random
import shutil
import subprocess
import sys
import time

from . import fasta, oracle

PROGDP = "--graphclustermethod none --graphtracemethod progdp"
# Variant names must not contain ":" (experiment.py splits "name:args" on the first colon).
VARIANTS = [
    ("default", ""),
    ("fm+opt", "--graphtracemethod fm --graphtraceoptimize true"),
    ("progdp", PROGDP + " --gcmx-order upgma"),
    # Soft constraints with MAGUS's normal evidence ("soft-m2", "soft-m3", "soft-m2-random" on
    # -b {bb}) and default+opt / progdp+opt were run on the first replicates and dropped: no gain.
    # richer evidence: the same backbones HMM-extended to all sequences (MAGUS "Slow" evidence)
    ("slow", "-b {ext}"),
    ("slow-soft-m2", "-s {rep}/split_m2 -b {ext}"),
    ("slow-soft-m3", "-s {rep}/split_m3 -b {ext}"),
    ("slow-soft-m4", "-s {rep}/split_m4 -b {ext}"),
    ("slow-soft-m3-random", "-s {rep}/split_m3_random -b {ext}"),
    # self-derived evidence: MAGUS's own output restricted to 8 random sequences per subset,
    # HMM-extended to all sequences (HMMER only, no new MAFFT runs)
    ("self-soft-m3", "-s {rep}/split_m3 -b {self}"),
]
ORACLES = [
    ("oracle-estSub-trueBB", "-b {o}/true_backbones"),
    ("oracle-estSub-trueFull", "-b {o}/true_full"),
    ("oracle-trueSub-estBB", "-s {o}/true_subalignments"),
    ("oracle-soft-m3-trueFull", "-s {rep}/split_m3 -b {o}/true_full"),
]
SPLITS = ((2, "linkage"), (3, "linkage"), (4, "linkage"), (2, "random"), (3, "random"))


def self_backbones(alignment, subalignment_dir, outdir, n=10, per_subset=8, seed=7):
    """Pseudo-backbones: the alignment restricted to `per_subset` random sequences of each subset."""
    aln = fasta.upper(fasta.read(alignment))
    subsets = [list(fasta.read(os.path.join(subalignment_dir, f))) for f in sorted(os.listdir(subalignment_dir))]
    rng = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)
    for b in range(n):
        taxa = [t for s in subsets for t in rng.sample(s, min(per_subset, len(s)))]
        fasta.write(fasta.restrict(aln, taxa), os.path.join(outdir, "self_backbone_{}.txt".format(b + 1)))


def main():
    results_dir, replicates = sys.argv[1], sys.argv[2:]
    os.makedirs(results_dir, exist_ok=True)
    results_path = os.path.join(results_dir, "results.jsonl")
    done = set()
    if os.path.exists(results_path):
        with open(results_path) as f:
            # only successful runs count as done, so failed variants (e.g. out of memory) are retried
            done = {(r["dataset"], r["variant"]) for r in map(json.loads, f) if "avgErr" in r}

    # free disk: drop working files (graphs are ~1 GB) of variants that already finished
    for log in glob.glob(os.path.join(results_dir, "*", "*", "log.txt")):
        workdir = os.path.dirname(log)
        if "finished in" in open(log, errors="replace").read():
            for item in os.listdir(workdir):
                if item not in ("log.txt", "stdout.log"):
                    path = os.path.join(workdir, item)
                    shutil.rmtree(path, ignore_errors=True) if os.path.isdir(path) else os.remove(path)

    for rep in replicates:
        rep = os.path.abspath(rep)
        dataset = os.path.basename(rep)
        inputs = os.path.join(rep, "inputs")
        oracle_dir = os.path.join(rep, "oracle")
        if not os.path.exists(os.path.join(oracle_dir, "true_full")):
            oracle.main(os.path.join(rep, "true.fasta"), os.path.join(inputs, "subalignments"),
                        os.path.join(inputs, "backbones"), oracle_dir)
        code_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        timing = {"dataset": dataset, "variant": "_timing"}  # seconds of the extra (non-merge) steps
        start = time.time()
        for m, method in SPLITS:
            split_dir = os.path.join(rep, "split_m{}{}".format(m, "_random" if method == "random" else ""))
            if not os.path.exists(split_dir):
                subprocess.run([sys.executable, "-m", "gcmx.split", os.path.join(inputs, "subalignments"),
                                split_dir, str(m), "--method", method], check=True, cwd=code_dir)
        timing["split_seconds_all_variants"] = round(time.time() - start, 1)
        ext_dir = os.path.join(rep, "ext_backbones")
        if not os.path.exists(ext_dir + ".DONE"):  # marker outside: MAGUS reads every file in a -b dir
            start = time.time()
            subprocess.run([sys.executable, "-m", "gcmx.extend", os.path.join(inputs, "backbones"),
                            os.path.join(rep, "unaligned.fasta"), ext_dir,
                            "--jobs", str(os.cpu_count())], check=True, cwd=code_dir)
            open(ext_dir + ".DONE", "w").close()
            timing["extend_seconds"] = round(time.time() - start, 1)
            timing["extend_jobs"] = os.cpu_count()
            with open(results_path, "a") as f:
                f.write(json.dumps(timing) + "\n")
        self_dir = os.path.join(rep, "ext_self")
        default_out = os.path.join(results_dir, dataset, "default.fasta")
        if not os.path.exists(self_dir + ".DONE"):
            if not os.path.exists(default_out):
                subprocess.run([sys.executable, "-m", "gcmx.experiment", "--dataset", dataset,
                                "--true", os.path.join(rep, "true.fasta"),
                                "--subalignments", os.path.join(inputs, "subalignments"),
                                "--backbones", os.path.join(inputs, "backbones"),
                                "--outdir", results_dir, "--variant", "default:"], check=True, cwd=code_dir)
                done.add((dataset, "default"))
            start = time.time()
            self_backbones(default_out, os.path.join(inputs, "subalignments"), os.path.join(rep, "self_bb"))
            subprocess.run([sys.executable, "-m", "gcmx.extend", os.path.join(rep, "self_bb"),
                            os.path.join(rep, "unaligned.fasta"), self_dir, "--jobs", str(os.cpu_count())],
                           check=True, cwd=code_dir)
            open(self_dir + ".DONE", "w").close()
            with open(results_path, "a") as f:
                f.write(json.dumps({"dataset": dataset, "variant": "_timing",
                                    "self_evidence_seconds": round(time.time() - start, 1)}) + "\n")
        fields = {"o": oracle_dir, "rep": rep, "bb": os.path.join(inputs, "backbones"), "ext": ext_dir,
                  "self": self_dir}
        only = set(filter(None, os.environ.get("PILOT_ONLY", "").split(",")))
        todo = [(name, extra.format(**fields)) for name, extra in VARIANTS + ORACLES
                if (dataset, name) not in done and (not only or name in only)]
        if not todo:
            continue
        cmd = [sys.executable, "-m", "gcmx.experiment", "--dataset", dataset,
               "--true", os.path.join(rep, "true.fasta"),
               "--subalignments", os.path.join(inputs, "subalignments"),
               "--backbones", os.path.join(inputs, "backbones"),
               "--outdir", results_dir, "--jobs", os.environ.get("PILOT_JOBS", "2")]
        for name, extra in todo:
            cmd += ["--variant", "{}:{}".format(name, extra)]
        subprocess.run(cmd, check=True, cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


if __name__ == "__main__":
    main()
