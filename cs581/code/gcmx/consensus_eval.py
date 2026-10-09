"""Evaluate the GCM consensus on the MAGUS paper's published alignments.

    python -m gcmx.consensus_eval OUT.jsonl WORKDIR DATASET/REP [DATASET/REP ...] [--jobs 3]
    (DATASET/REP as in Results.zip, e.g. 1000M2/R0, RNASim/R3, balibase/RV100_BBA0039)

Scenarios (inputs are the authors' published alignments of the same replicate):
  methods          MAGUS(Fast), MAGUS(Slow), PASTA(3), PASTA(3)+GCM
  pasta-iterations PASTA after 1, 3 and 4 iterations
For each replicate and scenario: score every input, the consensus, and record
which input the consensus chose as primary (the "most central" input), so the
consensus can be compared with the best input (oracle) and the central input
(what a user could pick without a reference). Resumable.
"""

import argparse
import concurrent.futures
import io
import json
import os
import shutil
import threading
import zipfile

from . import consensus, fasta, score
from .validate_published import URL, HttpFile

SCENARIOS = {
    "methods": ["gcm.txt", "gcm_slow.txt", "pasta_align.txt", "pasta_3_gcm_align.txt"],
    "pasta-iterations": ["pasta_1_align.txt", "pasta_align.txt", "pasta_4_align.txt"],
}
_zip_lock = threading.Lock()


def parse(text):
    seqs, name = {}, None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line and name is not None:
            seqs[name].append(line)
    return {n: "".join(s).upper() for n, s in seqs.items()}


def run_replicate(z, names, rep_path, workdir, done, out_path):
    key = rep_path.replace("/", "_")
    w = os.path.join(workdir, key)
    os.makedirs(w, exist_ok=True)
    files = {}
    with _zip_lock:
        for f in set(sum(SCENARIOS.values(), [])) | {"true_align.txt", "true_align_clean.txt"}:
            member = "Outputs/{}/{}".format(rep_path, f)
            if member in names:
                files[f] = parse(z.read(member).decode())
    est_taxa = set(files["gcm.txt"])
    ref = next(files[r] for r in ("true_align_clean.txt", "true_align.txt") if r in files and set(files[r]) == est_taxa)
    ref_path = os.path.join(w, "reference.fasta")
    fasta.write(ref, ref_path)
    paths = {}
    for f in files:
        if f.startswith("true"):
            continue
        paths[f] = os.path.join(w, f.replace(".txt", ".fasta"))
        fasta.write(files[f], paths[f])

    rows = []
    for scenario, inputs in SCENARIOS.items():
        if (rep_path, scenario) in done or not all(f in paths for f in inputs):
            continue
        row = {"replicate": rep_path, "scenario": scenario, "inputs": inputs}
        row["input_err"] = {f: score.fastsp(ref_path, paths[f])["avgErr"] for f in inputs}
        out = os.path.join(w, scenario + "_consensus.fasta")
        try:
            info = consensus.build(out, os.path.join(w, scenario), [paths[f] for f in inputs])
            row.update(info)
            row["central_input"] = inputs[info["primary"]]
            row["consensus"] = score.fastsp(ref_path, out)
        except Exception as e:  # recorded, not fatal
            row["error"] = repr(e)[:300]
        rows.append(row)
        with _zip_lock:
            with open(out_path, "a") as f:
                f.write(json.dumps(row) + "\n")
        print(json.dumps({k: row.get(k) for k in ("replicate", "scenario", "central_input")} |
                         {"best": min(row["input_err"].values()),
                          "consensus": row.get("consensus", {}).get("avgErr")}), flush=True)
    shutil.rmtree(os.path.join(w, "methods", "gcm"), ignore_errors=True)
    shutil.rmtree(os.path.join(w, "pasta-iterations", "gcm"), ignore_errors=True)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("out")
    parser.add_argument("workdir")
    parser.add_argument("replicates", nargs="+")
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    done = set()
    if os.path.exists(args.out):
        done = {(r["replicate"], r["scenario"]) for r in map(json.loads, open(args.out)) if "consensus" in r}
    z = zipfile.ZipFile(io.BufferedReader(HttpFile(URL), buffer_size=1 << 22))
    names = set(z.namelist())
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_replicate, z, names, r, args.workdir, done, args.out) for r in args.replicates]
        for fut in futures:
            fut.result()


if __name__ == "__main__":
    main()
