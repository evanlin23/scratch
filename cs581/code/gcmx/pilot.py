"""Pilot study: run all merge variants + oracle conditions on prepared replicates.

    python -m gcmx.pilot RESULTS_DIR REPLICATE_DIR [REPLICATE_DIR ...]

Each REPLICATE_DIR comes from gcmx.prep (true.fasta, inputs/subalignments,
inputs/backbones). (dataset, variant) pairs already in
RESULTS_DIR/results.jsonl are skipped, so the script can be rerun as more
replicates finish. PILOT_ONLY=name1,name2 restricts the variants (for large
datasets); PILOT_JOBS sets how many variants run concurrently.
"""

import json
import os
import subprocess
import sys

from . import oracle

PROGDP = "--graphclustermethod none --graphtracemethod progdp"
VARIANTS = [
    ("default", ""),
    ("default+opt", "--graphtraceoptimize true"),
    ("fm+opt", "--graphtracemethod fm --graphtraceoptimize true"),
    ("progdp", PROGDP + " --gcmx-order upgma"),
    ("progdp-grow", PROGDP + " --gcmx-order grow"),
    ("progdp+opt", PROGDP + " --gcmx-order upgma --graphtraceoptimize true"),
    ("weight-frac", "--gcmx-weight frac"),
    # soft constraints: split each subset alignment into m similarity groups, merge all groups at once
    ("soft-m2", "-s {rep}/split_m2 -b {bb}"),
    ("soft-m3", "-s {rep}/split_m3 -b {bb}"),
    ("soft-m2-random", "-s {rep}/split_m2_random -b {bb}"),  # control: random instead of similar groups
]
ORACLES = [
    ("oracle:estSub+trueBB", "-b {o}/true_backbones"),
    ("oracle:estSub+trueFull", "-b {o}/true_full"),
    ("oracle:trueSub+estBB", "-s {o}/true_subalignments"),
    ("oracle:trueSub+estBB+progdp", "-s {o}/true_subalignments " + PROGDP),
    ("oracle:soft-m2+trueFull", "-s {rep}/split_m2 -b {o}/true_full"),
]


def main():
    results_dir, replicates = sys.argv[1], sys.argv[2:]
    os.makedirs(results_dir, exist_ok=True)
    results_path = os.path.join(results_dir, "results.jsonl")
    done = set()
    if os.path.exists(results_path):
        with open(results_path) as f:
            done = {(r["dataset"], r["variant"]) for r in map(json.loads, f)}

    for rep in replicates:
        rep = os.path.abspath(rep)
        dataset = os.path.basename(rep)
        inputs = os.path.join(rep, "inputs")
        oracle_dir = os.path.join(rep, "oracle")
        if not os.path.exists(os.path.join(oracle_dir, "true_full")):
            oracle.main(os.path.join(rep, "true.fasta"), os.path.join(inputs, "subalignments"),
                        os.path.join(inputs, "backbones"), oracle_dir)
        for m, method in ((2, "linkage"), (3, "linkage"), (2, "random")):
            split_dir = os.path.join(rep, "split_m{}{}".format(m, "_random" if method == "random" else ""))
            if not os.path.exists(split_dir):
                subprocess.run([sys.executable, "-m", "gcmx.split", os.path.join(inputs, "subalignments"),
                                split_dir, str(m), "--method", method], check=True,
                               cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        fields = {"o": oracle_dir, "rep": rep, "bb": os.path.join(inputs, "backbones")}
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
