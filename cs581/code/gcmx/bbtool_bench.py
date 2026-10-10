"""Does the backbone aligner matter? MAGUS with its GCM backbones aligned by another tool.

    python -m gcmx.bbtool_bench JOBFILE OUT.jsonl WORKDIR [--draws 0,1,2] [--tools clustalo,mafft-auto]
                                [--e2e clustalo] [--threads 4]

JOBFILE lines: NAME TRUE_ALIGNMENT NUM_SUBSETS (as in fanout/jobs_*.txt). MAGUS is unseeded, so
every draw is a fresh MAGUS run (new decomposition and new backbone sequence sets). For every
draw, then every job, sequentially, all with the same thread count:

  magus         MAGUS(Fast) end to end from unaligned sequences, the paper's flags (timed)
  e2e-X         MAGUS end to end with its backbones aligned by X (gcmx.run_magus --gcmx-backbonetool X;
                timed; its own decomposition and backbone draw, so unpaired with `magus`)
  merge-mafft   control: GCM merge only on magus's own subsets and backbones (must reproduce `magus`)
  merge-X       paired: magus's own subsets and the same backbone sequence sets, realigned by X,
                then the GCM merge only. Records backbone alignment time, merge time, the
                backbones' own SPFN/SPFP against the reference, and how many residue pairs they align.

X is a key of run_magus.BACKBONE_TOOLS (clustalo, mafft-auto, linsi-noep, ginsi). One row per
(dataset, draw) is appended to OUT.jsonl when all requested methods are done; partial progress is
kept in WORKDIR/<name>_d<draw>/state.json, so a killed run resumes where it stopped (rerun the same
command; add tools later the same way).
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time

from . import fasta, score
from .e2e_bench import CODE, REPO, acc, is_protein, magus_flags, timed
from .run_magus import BACKBONE_TOOLS

MERGE_FLAGS = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
               "-f", "4"]


def mafft_path():
    from magus.configuration import Configs
    return Configs.mafftPath


def align_backbones(tool, src_dir, dst_dir, threads, log_dir):
    """Realign every backbone_N_unalign.txt in src_dir with `tool`, one after another with `threads`
    threads each (what MAGUS does with -np threads). Returns wall-clock seconds."""
    os.makedirs(dst_dir, exist_ok=True)
    start = time.time()
    for f in sorted(os.listdir(src_dir)):
        if not f.endswith("_unalign.txt"):
            continue
        argv = BACKBONE_TOOLS[tool](threads)
        if argv[0] == "mafft":
            argv[0] = mafft_path()
        out = os.path.join(dst_dir, f.replace("_unalign.txt", "_mafft.txt"))
        with open(out, "w") as o, open(os.path.join(log_dir, "bb_{}_{}.log".format(tool, f)), "w") as e:
            subprocess.run(argv + [os.path.join(src_dir, f)], stdout=o, stderr=e, check=True)
    return round(time.time() - start, 1)


def backbone_stats(ref, bb_dir):
    """Mean SPFN/SPFP of the backbone alignments against the reference restricted to their
    sequences, and the mean number of aligned residue pairs (= evidence pairs GCM receives)."""
    fn, fp, pairs, n = 0.0, 0.0, 0, 0
    for f in sorted(os.listdir(bb_dir)):
        if not f.endswith("_mafft.txt"):
            continue
        est = fasta.upper(fasta.read(os.path.join(bb_dir, f)))
        tr = os.path.join(bb_dir, "true_" + f)
        fasta.write(fasta.restrict(ref, list(est)), tr)
        s = score.fastsp(tr, os.path.join(bb_dir, f))
        os.remove(tr)
        fn, fp, n = fn + s["SPFN"], fp + s["SPFP"], n + 1
        cols = zip(*est.values())
        pairs += sum(k * (k - 1) // 2 for k in (sum(c != "-" for c in col) for col in cols))
    return {"bb_SPFN": round(fn / n, 4), "bb_SPFP": round(fp / n, 4), "bb_pairs": pairs // n}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("jobs")
    parser.add_argument("out")
    parser.add_argument("work")
    parser.add_argument("--draws", default="0,1,2")
    parser.add_argument("--tools", default="clustalo,mafft-auto", help="paired merge-only backbone tools")
    parser.add_argument("--e2e", default="clustalo", help="end-to-end MAGUS with these backbone tools ('' = none)")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    tools = [t for t in args.tools.split(",") if t]
    e2e = [t for t in args.e2e.split(",") if t]
    for t in tools + e2e:
        if t not in BACKBONE_TOOLS or t == "mafft":
            raise SystemExit("unknown backbone tool: " + t)
    done = set()
    if os.path.exists(args.out):
        done = {(r["dataset"], r["draw"]) for r in map(json.loads, open(args.out))}
    T = str(args.threads)
    py = [sys.executable, "-m"]
    jobs = [l.split()[:3] for l in open(args.jobs) if l.strip()]

    for draw in map(int, args.draws.split(",")):
        for name, src, k in jobs:
            if (name, draw) in done:
                continue
            src = src if os.path.isabs(src) else os.path.join(REPO, src)
            w = os.path.join(args.work, "{}_d{}".format(name, draw))
            os.makedirs(w, exist_ok=True)
            state_path = os.path.join(w, "state.json")
            row = json.load(open(state_path)) if os.path.exists(state_path) else {}
            ref = fasta.upper(fasta.read(src))
            true, unaligned = os.path.join(w, "true.fasta"), os.path.join(w, "unaligned.fasta")
            fasta.write(ref, true)
            fasta.write(fasta.ungap(ref), unaligned)
            row.update({"dataset": name, "draw": draw, "threads": args.threads, "nproc": os.cpu_count(),
                        "nseq": len(ref), "protein": is_protein(ref)})

            def log_row(method, data):
                row[method] = data
                with open(state_path + ".tmp", "w") as f:
                    json.dump(row, f)
                os.replace(state_path + ".tmp", state_path)
                print(json.dumps({name: {"draw": draw, method: data}}), flush=True)

            # baseline MAGUS; its subsets and backbone sequence sets are kept in w/inputs for the paired merges
            inputs = os.path.join(w, "inputs")
            if "magus" not in row or not os.path.isdir(inputs):
                shutil.rmtree(os.path.join(w, "magus"), ignore_errors=True)
                shutil.rmtree(inputs, ignore_errors=True)
                out = os.path.join(w, "magus.fasta")
                wall, cpu = timed(py + ["gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", T, "-d",
                                        os.path.join(w, "magus"), "-i", unaligned, "-o", out] + magus_flags(k),
                                  os.path.join(w, "magus.log"))
                shutil.copytree(os.path.join(w, "magus", "subalignments"), os.path.join(inputs, "subalignments"))
                os.makedirs(os.path.join(inputs, "backbones"))
                for f in os.listdir(os.path.join(w, "magus", "graph")):
                    if f.startswith("backbone_") and f.endswith(("_unalign.txt", "_mafft.txt")):
                        shutil.copy(os.path.join(w, "magus", "graph", f), os.path.join(inputs, "backbones"))
                shutil.rmtree(os.path.join(w, "magus"), ignore_errors=True)
                log_row("magus", {"wall": wall, "cpu": cpu, **acc(true, out),
                                  **backbone_stats(ref, os.path.join(inputs, "backbones"))})

            for tool in e2e:
                key = "e2e-" + tool
                if key in row:
                    continue
                d = os.path.join(w, key)
                shutil.rmtree(d, ignore_errors=True)
                out = os.path.join(w, key + ".fasta")
                wall, cpu = timed(py + ["gcmx.run_magus", "--gcmx-fastgraph", "false", "--gcmx-backbonetool", tool,
                                        "-np", T, "-d", d, "-i", unaligned, "-o", out] + magus_flags(k),
                                  os.path.join(w, key + ".log"))
                stats = backbone_stats(ref, os.path.join(d, "graph"))
                shutil.rmtree(d, ignore_errors=True)
                log_row(key, {"wall": wall, "cpu": cpu, **acc(true, out), **stats})

            for tool in ["mafft"] + tools:
                key = "merge-" + tool
                if key in row:
                    continue
                d = os.path.join(w, key)
                shutil.rmtree(d, ignore_errors=True)
                os.makedirs(d)
                if tool == "mafft":
                    bb, bb_wall = os.path.join(inputs, "backbones"), None
                else:
                    bb = os.path.join(d, "backbones")
                    bb_wall = align_backbones(tool, os.path.join(inputs, "backbones"), bb, args.threads, d)
                bb_only = os.path.join(d, "bb_aligned")
                os.makedirs(bb_only)
                for f in os.listdir(bb):
                    if f.endswith("_mafft.txt"):
                        shutil.copy(os.path.join(bb, f), bb_only)
                out = os.path.join(w, key + ".fasta")
                m_wall, m_cpu = timed(py + ["gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", T, "-d",
                                            os.path.join(d, "magus"), "-s", os.path.join(inputs, "subalignments"),
                                            "-b", bb_only, "-o", out] + MERGE_FLAGS, os.path.join(w, key + ".log"))
                data = {"merge_wall": m_wall, "merge_cpu": m_cpu, **acc(true, out), **backbone_stats(ref, bb_only)}
                if bb_wall is not None:
                    data["backbone_wall"] = bb_wall
                shutil.rmtree(d, ignore_errors=True)
                log_row(key, data)

            with open(args.out, "a") as f:
                f.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()
