"""Does the backbone aligner matter? MAGUS with its GCM backbones aligned by another tool.

    python -m gcmx.bbtool_bench JOBFILE OUT.jsonl WORKDIR [--draws 0,1,2] [--tools clustalo,mafft-auto]
                                [--e2e clustalo] [--threads 4]

JOBFILE lines: NAME TRUE_ALIGNMENT NUM_SUBSETS [UNALIGNED] (as in fanout/jobs_*.txt). With the optional
UNALIGNED file (e.g. HomFam: the whole family, while the reference covers only the Homstrad seeds), MAGUS
aligns UNALIGNED and every score is computed on the estimated alignment restricted to the reference's
sequences (all-gap columns removed); backbone stats use only the reference sequences each backbone holds. MAGUS is unseeded, so
every draw is a fresh MAGUS run (new decomposition and new backbone sequence sets). For every
draw, then every job, sequentially, all with the same thread count:

  magus         MAGUS(Fast) end to end from unaligned sequences, the paper's flags (timed)
  e2e-X         MAGUS end to end with its backbones aligned by X (gcmx.run_magus --gcmx-backbonetool X;
                timed; its own decomposition and backbone draw, so unpaired with `magus`)
  merge-mafft   control: GCM merge only on magus's own subsets and backbones (must reproduce `magus`)
  merge-X       paired: magus's own subsets and the same backbone sequence sets, realigned by X,
                then the GCM merge only. Records backbone alignment time, merge time, the
                backbones' own SPFN/SPFP against the reference, and how many residue pairs they align.
  merge-union-X paired: magus's 10 L-INS-i backbones plus the same 10 sequence sets realigned by X
                (20 backbones, two aligners: does mixing error-decorrelated evidence help?)

X is a key of run_magus.BACKBONE_TOOLS (clustalo, mafft-auto, linsi-noep, ginsi). One row per
(dataset, draw) is appended to OUT.jsonl when all requested methods are done; partial progress is
kept in WORKDIR/<name>_d<draw>/state.json, so a killed run resumes where it stopped (rerun the same
command; add tools later the same way).
"""

import argparse
import concurrent.futures
import json
import os
import shutil
import subprocess
import sys
import time

from . import fasta, score
from .e2e_bench import CODE, REPO, is_protein, magus_flags, timed
from .run_magus import BACKBONE_TOOLS

# Bumped when timings of e2e-/merge- rows change meaning; older rows are redone on the next run
# (v2: Clustal Omega single-threaded and backbones realigned in parallel, as MAGUS schedules them).
VERSION = 2
# GCMX_FASTGRAPH=true: the vectorized graph builder (gcmx.fastgraph; same graph, faster)
FASTGRAPH = os.environ.get("GCMX_FASTGRAPH", "false")
MERGE_FLAGS = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
               "-f", "4"]


def mafft_path():
    from magus.configuration import Configs
    return Configs.mafftPath


def align_backbones(tool, src_dir, dst_dir, threads, log_dir):
    """Realign every backbone_N_unalign.txt in src_dir with `tool` the way MAGUS schedules its
    backbone tasks with -np `threads`: up to `threads` at once, each given BACKBONE_TOOLS[tool](threads)
    (MAFFT gets --thread threads, Clustal Omega 1 thread). Returns wall-clock seconds."""
    os.makedirs(dst_dir, exist_ok=True)

    def one(f):
        argv = BACKBONE_TOOLS[tool](threads)
        if argv[0] == "mafft":
            argv[0] = mafft_path()
        out = os.path.join(dst_dir, f.replace("_unalign.txt", "_mafft.txt"))
        with open(out, "w") as o, open(os.path.join(log_dir, "bb_{}_{}.log".format(tool, f)), "w") as e:
            subprocess.run(argv + [os.path.join(src_dir, f)], stdout=o, stderr=e, check=True)

    start = time.time()
    files = [f for f in sorted(os.listdir(src_dir)) if f.endswith("_unalign.txt")]
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        list(pool.map(one, files))
    return round(time.time() - start, 1)


def acc_ref(true, path):
    """SPFN/SPFP/TC of `path` against `true`, on the estimate restricted to the reference's sequences
    (a no-op when both hold the same sequences)."""
    ref = fasta.read(true)
    est = fasta.read(path)
    if set(est) != set(ref):
        tmp = path + ".refonly.fasta"
        fasta.write(fasta.restrict(est, list(ref)), tmp)
        path = tmp
    s = score.fastsp(true, path)
    return {k: s[k] for k in ("SPFN", "SPFP", "avgErr", "TC", "LenEst", "LenRef") if k in s}


def backbone_stats(ref, bb_dir):
    """Mean SPFN/SPFP of the backbone alignments against the reference restricted to their
    sequences, and the mean number of aligned residue pairs (= evidence pairs GCM receives)."""
    fn, fp, pairs, n = 0.0, 0.0, 0, 0
    for f in sorted(os.listdir(bb_dir)):
        if not f.endswith("_mafft.txt"):
            continue
        est = fasta.upper(fasta.read(os.path.join(bb_dir, f)))
        cols = zip(*est.values())
        pairs += sum(k * (k - 1) // 2 for k in (sum(c != "-" for c in col) for col in cols))
        shared = [t for t in est if t in ref]
        if len(shared) < 2:
            continue
        tr, es = os.path.join(bb_dir, "true_" + f), os.path.join(bb_dir, "est_" + f)
        fasta.write(fasta.restrict(ref, shared), tr)
        fasta.write(fasta.restrict(est, shared), es)
        s = score.fastsp(tr, es)
        os.remove(tr)
        os.remove(es)
        fn, fp, n = fn + s["SPFN"], fp + s["SPFP"], n + 1
    nb = len([f for f in os.listdir(bb_dir) if f.endswith("_mafft.txt")])
    if n == 0:
        return {"bb_SPFN": None, "bb_SPFP": None, "bb_pairs": pairs // max(nb, 1), "bb_scored": 0}
    return {"bb_SPFN": round(fn / n, 4), "bb_SPFP": round(fp / n, 4), "bb_pairs": pairs // nb, "bb_scored": n}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("jobs")
    parser.add_argument("out")
    parser.add_argument("work")
    parser.add_argument("--draws", default="0,1,2")
    parser.add_argument("--tools", default="clustalo,mafft-auto", help="paired merge-only backbone tools")
    parser.add_argument("--e2e", default="clustalo", help="end-to-end MAGUS with these backbone tools ('' = none)")
    parser.add_argument("--union", default="clustalo", help="merge-union-X for these X in --tools ('' = none)")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    tools = [t for t in args.tools.split(",") if t]
    unions = [t for t in args.union.split(",") if t in tools]
    e2e = [t for t in args.e2e.split(",") if t]
    for t in tools + e2e:
        if t not in BACKBONE_TOOLS or t == "mafft":
            raise SystemExit("unknown backbone tool: " + t)
    done = set()
    if os.path.exists(args.out):
        done = {(r["dataset"], r["draw"]) for r in map(json.loads, open(args.out))}
    T = str(args.threads)
    py = [sys.executable, "-m"]
    jobs = [(l.split() + [None])[:4] for l in open(args.jobs) if l.strip()]

    for draw in map(int, args.draws.split(",")):
        for name, src, k, unaln_src in jobs:
            if (name, draw) in done:
                continue
            src = src if os.path.isabs(src) else os.path.join(REPO, src)
            w = os.path.join(args.work, "{}_d{}".format(name, draw))
            os.makedirs(w, exist_ok=True)
            state_path = os.path.join(w, "state.json")
            row = json.load(open(state_path)) if os.path.exists(state_path) else {}
            ref = fasta.upper(fasta.read(src))
            true, unaligned = os.path.join(w, "true.fasta"), os.path.join(w, "unaligned.fasta")
            if unaln_src:
                unaln_src = unaln_src if os.path.isabs(unaln_src) else os.path.join(REPO, unaln_src)
                seqs = fasta.upper(fasta.ungap(fasta.read(unaln_src)))
                missing = [t for t in ref if t not in seqs]
                if missing:
                    raise SystemExit("{}: {} reference sequences not in {}".format(name, len(missing), unaln_src))
                fasta.write(ref, true)
                fasta.write(seqs, unaligned)
            else:
                seqs = fasta.ungap(ref)
                fasta.write(ref, true)
                fasta.write(seqs, unaligned)
            row.update({"dataset": name, "draw": draw, "threads": args.threads, "nproc": os.cpu_count(),
                        "nseq": len(seqs), "nref": len(ref), "protein": is_protein(ref)})

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
                if os.path.exists(out):
                    os.remove(out)
                wall, cpu = timed(py + ["gcmx.run_magus", "--gcmx-fastgraph", FASTGRAPH, "-np", T, "-d",
                                        os.path.join(w, "magus"), "-i", unaligned, "-o", out] + magus_flags(k),
                                  os.path.join(w, "magus.log"))
                shutil.copytree(os.path.join(w, "magus", "subalignments"), os.path.join(inputs, "subalignments"))
                os.makedirs(os.path.join(inputs, "backbones"))
                for f in os.listdir(os.path.join(w, "magus", "graph")):
                    if f.startswith("backbone_") and f.endswith(("_unalign.txt", "_mafft.txt")):
                        shutil.copy(os.path.join(w, "magus", "graph", f), os.path.join(inputs, "backbones"))
                shutil.rmtree(os.path.join(w, "magus"), ignore_errors=True)
                log_row("magus", {"wall": wall, "cpu": cpu, **acc_ref(true, out),
                                  **backbone_stats(ref, os.path.join(inputs, "backbones"))})

            for tool in e2e:
                key = "e2e-" + tool
                if row.get(key, {}).get("v", 1) >= VERSION:
                    continue
                d = os.path.join(w, key)
                shutil.rmtree(d, ignore_errors=True)
                out = os.path.join(w, key + ".fasta")
                if os.path.exists(out):
                    os.remove(out)  # MAGUS skips the whole run when its output file exists
                wall, cpu = timed(py + ["gcmx.run_magus", "--gcmx-fastgraph", FASTGRAPH, "--gcmx-backbonetool", tool,
                                        "-np", T, "-d", d, "-i", unaligned, "-o", out] + magus_flags(k),
                                  os.path.join(w, key + ".log"))
                stats = backbone_stats(ref, os.path.join(d, "graph"))
                shutil.rmtree(d, ignore_errors=True)
                log_row(key, {"wall": wall, "cpu": cpu, **acc_ref(true, out), **stats, "v": VERSION})

            def realigned(tool, d):
                """X-realigned backbones, kept in inputs/bb_X (small) for merge-union-X; returns (dir, wall)."""
                bb = os.path.join(inputs, "bb_" + tool)
                if os.path.exists(os.path.join(bb, "done.json")):
                    return bb, json.load(open(os.path.join(bb, "done.json")))["wall"]
                shutil.rmtree(bb, ignore_errors=True)
                wall = align_backbones(tool, os.path.join(inputs, "backbones"), bb, args.threads, d)
                json.dump({"wall": wall}, open(os.path.join(bb, "done.json"), "w"))
                return bb, wall

            for tool in ["mafft"] + tools + ["union-" + t for t in unions]:
                key = "merge-" + tool
                if row.get(key, {}).get("v", 1) >= VERSION:
                    continue
                d = os.path.join(w, key)
                shutil.rmtree(d, ignore_errors=True)
                os.makedirs(d)
                bb_only = os.path.join(d, "bb_aligned")
                os.makedirs(bb_only)
                if tool == "mafft":
                    sources, bb_wall = [("", os.path.join(inputs, "backbones"))], None
                elif tool.startswith("union-"):
                    bb, bb_wall = realigned(tool[len("union-"):], d)
                    sources = [("linsi_", os.path.join(inputs, "backbones")), (tool[len("union-"):] + "_", bb)]
                else:
                    bb, bb_wall = realigned(tool, d)
                    sources = [("", bb)]
                for prefix, src in sources:
                    for f in os.listdir(src):
                        if f.endswith("_mafft.txt"):
                            shutil.copy(os.path.join(src, f), os.path.join(bb_only, prefix + f))
                out = os.path.join(w, key + ".fasta")
                if os.path.exists(out):
                    os.remove(out)
                m_wall, m_cpu = timed(py + ["gcmx.run_magus", "--gcmx-fastgraph", FASTGRAPH, "-np", T, "-d",
                                            os.path.join(d, "magus"), "-s", os.path.join(inputs, "subalignments"),
                                            "-b", bb_only, "-o", out] + MERGE_FLAGS, os.path.join(w, key + ".log"))
                data = {"merge_wall": m_wall, "merge_cpu": m_cpu, **acc_ref(true, out), **backbone_stats(ref, bb_only),
                        "v": VERSION}
                if bb_wall is not None:
                    data["backbone_wall"] = bb_wall
                shutil.rmtree(d, ignore_errors=True)
                log_row(key, data)

            with open(args.out, "a") as f:
                f.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()
