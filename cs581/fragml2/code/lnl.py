"""Re-evaluate the log-likelihood of final trees under one fixed model per replicate.

    python lnl.py [--datasets M1HF RNASimHF] [--arms ...] [--workers 4]

Model = that replicate's base_raxmlng GTR+G parameters (trees/<aln>.base_raxmlng.model), RAxML-NG
--evaluate --opt-model off (branch lengths optimised), 1 thread. Rows -> results/lnl.jsonl (skips done).
"""
import argparse
import json
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
RX = "/opt/mm/root/envs/fml/bin/raxml-ng"
OUT = os.path.join(HERE, "..", "results", "lnl.jsonl")


def ev(row):
    d = os.path.dirname(os.path.dirname(row["tree"]))
    model = os.path.join(d, "trees", row["aln"] + ".base_raxmlng.model")
    msa = os.path.join(d, row["aln"] + ".clean.fasta")
    w = tempfile.mkdtemp(prefix="lnl_")
    subprocess.run([RX, "--evaluate", "--msa", msa, "--tree", row["tree"], "--model", model, "--opt-model", "off",
                    "--threads", "1", "--prefix", os.path.join(w, "e"), "--redo"], stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, check=True)
    lnl = None
    for line in open(os.path.join(w, "e.raxml.log")):
        if line.startswith("Final LogLikelihood:"):
            lnl = float(line.split(":")[1])
    r = {k: row[k] for k in ("dataset", "rep", "aln", "arm")}
    r["lnl"] = lnl
    with open(OUT, "a") as f:
        f.write(json.dumps(r) + "\n")
    print(json.dumps(r), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["M1HF", "RNASimHF"])
    ap.add_argument("--arms", nargs="*")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    done = set()
    if os.path.exists(OUT):
        done = {(r["dataset"], r["rep"], r["aln"], r["arm"]) for r in map(json.loads, open(OUT))}
    rows = [json.loads(l) for l in open(os.path.join(HERE, "..", "results", "runs.jsonl"))]
    todo = []
    for r in rows:
        k = (r["dataset"], r["rep"], r["aln"], r["arm"])
        if r["dataset"] not in a.datasets or k in done or (a.arms and r["arm"] not in a.arms):
            continue
        d = os.path.dirname(os.path.dirname(r["tree"]))
        if os.path.exists(os.path.join(d, "trees", r["aln"] + ".base_raxmlng.model")):
            todo.append(r)
            done.add(k)
    with ThreadPoolExecutor(a.workers) as ex:
        list(ex.map(ev, todo))


if __name__ == "__main__":
    main()
