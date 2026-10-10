"""Fast MAGUS: measured end-to-end runtime + accuracy of cheaper guide trees / fewer backbones
(+ self-derived soft evidence) against MAGUS with the paper's flags, on one idle machine.

    python fm_bench.py JOBFILE OUT.jsonl WORKDIR [--threads 4]

JOBFILE lines:  NAME TRUE_ALIGNMENT NUM_SUBSETS VARIANT[,VARIANT...]
VARIANT names (see BASES): a base MAGUS run, optionally followed by "+ss" / "+ss2"
(one or two rounds of self-soft on that base run's output). Every base run starts from
the unaligned sequences; "+ss" rows add the measured wall-clock of their extra steps to
the base run's measured wall-clock (they literally run afterwards on its files).

Restartable per base run: state is kept in WORKDIR/NAME/state.json.
"""

import argparse
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CODE = os.path.join(REPO, "cs581", "code")  # gcmx + MAGUS live here (not modified)
sys.path.insert(0, CODE)

from gcmx import fasta, score  # noqa: E402
from gcmx.pilot import self_backbones  # noqa: E402

# name -> (guide tree type, skeleton size, backbones, "free" speedups (vectorized graph build,
# identical graph; MCL with 4 threads, identical clusters))
BASES = {
    "magus": ("fasttree", 300, 10, False),      # MAGUS(Fast), the paper's flags, pure MAGUS
    "magus-free": ("fasttree", 300, 10, True),
    "r4": ("fasttree", 300, 4, True),
    "r8": ("fasttree", 300, 8, True),
    "sk100-r4": ("fasttree", 100, 4, True),
    "sk100-r2": ("fasttree", 100, 2, True),
    "sk150-r4": ("fasttree", 150, 4, True),
    "pt-r4": ("parttree", 300, 4, True),
    "pt-r8": ("parttree", 300, 8, True),
    "pt-r10": ("parttree", 300, 10, True),
    "cl-r4": ("clustal", 300, 4, True),
    "cl-r10": ("clustal", 300, 10, True),
    "pt-r2": ("parttree", 300, 2, True),
}


def magus_flags(k, guide, skel, r):
    return ["--maxsubsetsize", "0", "--maxnumsubsets", str(k), "--decompstrategy", "pastastyle",
            "--decompskeletonsize", str(skel), "-t", guide, "--graphbuildmethod", "mafft",
            "--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
            "-r", str(r), "-m", "200", "-f", "4"]


def timed(cmd, log):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.time()
    with open(log, "w") as f:
        proc = subprocess.run(cmd, cwd=CODE, stdout=f, stderr=subprocess.STDOUT)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    if proc.returncode != 0:
        raise RuntimeError("failed ({}): {}".format(proc.returncode, " ".join(cmd)))
    return round(time.time() - start, 1), round(cpu, 1)


def acc(true, path):
    s = score.fastsp(true, path)
    return {k: s[k] for k in ("SPFN", "SPFP", "avgErr", "LenEst", "LenRef")}


STAMP = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)")
SECS = (("tree", r"Built initial tree on .* in ([0-9.e+-]+) sec"), ("decomp", r"Decomposed .* in ([0-9.e+-]+) sec"),
        ("graph", r"Built the alignment graph in ([0-9.e+-]+) sec"), ("cluster", r"Clustered the graph in ([0-9.e+-]+) sec"),
        ("trace", r"Found alignment graph trace in ([0-9.e+-]+) sec"), ("merge", r"Merged .* in ([0-9.e+-]+) sec"),
        ("total", r"MAGUS finished in ([0-9.e+-]+) sec"))


def stages(logpath):
    """Seconds per MAGUS stage from its log.txt, plus when the last subset alignment and the last
    backbone finished (seconds after start)."""
    out, first = {}, None
    if not os.path.exists(logpath):
        return out
    for line in open(logpath, errors="replace"):
        m = STAMP.match(line)
        if not m:
            continue
        t = time.mktime(time.strptime(m.group(1), "%Y-%m-%d %H:%M:%S"))
        first = t if first is None else first
        for k, pat in SECS:
            mm = re.search(pat, line)
            if mm:
                out[k] = round(float(mm.group(1)), 1)
        if "Completed a task" in line:
            for k, pat in (("subsets_done", "/subalignment_"), ("backbones_done", "/backbone_")):
                if pat in line:
                    out[k] = t - first
    return out


def clean(workdir):
    """Drop a finished MAGUS working dir (graphs are large); keep its log."""
    for item in os.listdir(workdir) if os.path.isdir(workdir) else []:
        if item != "log.txt":
            path = os.path.join(workdir, item)
            shutil.rmtree(path, ignore_errors=True) if os.path.isdir(path) else os.remove(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("jobs")
    parser.add_argument("out")
    parser.add_argument("work")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--groups", type=int, default=3)
    args = parser.parse_args()
    T = str(args.threads)
    py = [sys.executable, "-m"]

    for line in open(args.jobs):
        if not line.strip() or line.startswith("#"):
            continue
        name, src, k, variants = line.split()[:4]
        variants = variants.split(",")
        src = src if os.path.isabs(src) else os.path.join(REPO, src)
        w = os.path.join(args.work, name)
        os.makedirs(w, exist_ok=True)
        state_path = os.path.join(w, "state.json")
        row = json.load(open(state_path)) if os.path.exists(state_path) else {}
        if all(v in row for v in variants):
            continue
        ref = fasta.upper(fasta.read(src))
        true, unaligned = os.path.join(w, "true.fasta"), os.path.join(w, "unaligned.fasta")
        fasta.write(ref, true)
        fasta.write(fasta.ungap(ref), unaligned)
        row.update({"dataset": name, "threads": args.threads, "nproc": os.cpu_count(), "nseq": len(ref)})

        def save():
            with open(state_path + ".tmp", "w") as f:
                json.dump(row, f)
            os.replace(state_path + ".tmp", state_path)

        def log_row(method, data):
            row[method] = data
            save()
            print(json.dumps({name: {method: data}}), flush=True)

        bases = []
        for v in variants:
            b = v.split("+")[0]
            if b not in bases:
                bases.append(b)
        for b in bases:
            guide, skel, r, free = BASES[b]
            wanted = [v for v in variants if v.split("+")[0] == b]
            if all(v in row for v in wanted) or "error" in row.get(b, {}):
                continue
            bw = os.path.join(w, b)
            out = os.path.join(w, b + ".fasta")
            rounds = max([2 if v.endswith("+ss2") else 1 if v.endswith("+ss") else 0 for v in wanted])
            subs = os.path.join(bw, "subalignments")
            if b not in row or (rounds > 0 and not os.path.isdir(subs)):
                shutil.rmtree(bw, ignore_errors=True)
                for v in wanted:
                    row.pop(v, None)
                extra = ["--gcmx-fastgraph", "false"] if not free else ["--gcmx-mclthreads", T]
                try:
                    wall, cpu = timed(py + ["gcmx.run_magus"] + extra + ["-np", T, "-d", bw, "-i", unaligned,
                                                                        "-o", out] + magus_flags(k, guide, skel, r),
                                      os.path.join(w, b + ".log"))
                except RuntimeError as e:
                    log_row(b, {"error": str(e)})
                    continue
                log_row(b, {"wall": wall, "cpu": cpu, "stages": stages(os.path.join(bw, "log.txt")),
                            **acc(true, out)})
            if rounds == 0:
                clean(bw)
                continue
            split = os.path.join(w, b + "_split")
            shutil.rmtree(split, ignore_errors=True)
            s_wall, s_cpu = timed(py + ["gcmx.split", subs, split, str(args.groups)], os.path.join(w, b + "_split.log"))
            wall, cpu = row[b]["wall"] + s_wall, row[b]["cpu"] + s_cpu
            prev = out
            for rnd in range(1, rounds + 1):
                tag = b + ("+ss" if rnd == 1 else "+ss{}".format(rnd))
                bb, ev, mw = (os.path.join(w, tag + s) for s in ("_bb", "_ev", "_merge"))
                for d in (bb, ev, mw):
                    shutil.rmtree(d, ignore_errors=True)
                start = time.time()
                self_backbones(prev, subs, bb, seed=7 + rnd - 1)
                prep = time.time() - start
                e_wall, e_cpu = timed(py + ["gcmx.extend", bb, unaligned, ev, "--jobs", T],
                                      os.path.join(w, tag + "_extend.log"))
                rout = os.path.join(w, tag + ".fasta")
                m_wall, m_cpu = timed(py + ["gcmx.run_magus", "-np", T, "--gcmx-mclthreads", T, "-d", mw,
                                            "-s", split, "-b", ev, "-o", rout], os.path.join(w, tag + "_merge.log"))
                wall += e_wall + prep + m_wall
                cpu += e_cpu + m_cpu
                if tag in variants:
                    log_row(tag, {"wall": round(wall, 1), "cpu": round(cpu, 1), "split_wall": s_wall,
                                  "evidence_wall": round(e_wall + prep, 1), "merge_wall": m_wall, **acc(true, rout)})
                prev = rout
                for d in (bb, ev, mw):
                    shutil.rmtree(d, ignore_errors=True)
            shutil.rmtree(split, ignore_errors=True)
            clean(bw)

        with open(args.out, "a") as f:
            f.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()
