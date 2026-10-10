"""Full-dataset benchmark: run aligners on unaligned sequences, score with FastSP.

    python3 bench.py JOBS OUT.jsonl --tools a,b,... [--work /opt/runs/sota] [--timeout 3600]

JOBS lines: NAME TRUE_ALIGNMENT [K]   (K = MAGUS max subsets; default 25)
Tool "magus" runs MAGUS(Fast) with the paper's flags (gcmx.e2e_bench.magus_flags).
Restartable: (dataset, tool) pairs already in OUT are skipped. Jobs run one at a
time with 4 threads; run nothing else heavy alongside.
"""

import argparse
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
GCMX_PARENT = os.path.join(REPO, "cs581", "code")
sys.path.insert(0, GCMX_PARENT)
sys.path.insert(0, HERE)

from gcmx import fasta, score  # noqa: E402
import tools  # noqa: E402


def magus_tool(k):
    from gcmx.e2e_bench import magus_flags

    def build(i, o, t, w):
        d = os.path.join(w, "magus_work")
        shutil.rmtree(d, ignore_errors=True)
        return ([sys.executable, "-m", "gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", str(t), "-d", d,
                 "-i", i, "-o", o, "--graphbuildhmmextend", "false"] + magus_flags(k), None)
    return build


def main():
    p = argparse.ArgumentParser()
    p.add_argument("jobs")
    p.add_argument("out")
    p.add_argument("--tools", required=True)
    p.add_argument("--work", default="/opt/runs/sota")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--timeout", type=float, default=3600)
    p.add_argument("--keep", action="store_true", help="keep estimated alignments")
    a = p.parse_args()

    done = set()
    if os.path.exists(a.out):
        done = {(r["dataset"], r["tool"]) for r in map(json.loads, open(a.out))}
    for line in open(a.jobs):
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        name, src = parts[0], parts[1]
        k = int(parts[2]) if len(parts) > 2 else 25
        src = src if os.path.isabs(src) else os.path.join(REPO, src)
        w = os.path.join(a.work, name)
        os.makedirs(w, exist_ok=True)
        true, unal = os.path.join(w, "true.fasta"), os.path.join(w, "unaligned.fasta")
        if not os.path.exists(unal):
            ref = fasta.upper(fasta.read(src))
            fasta.write(ref, true)
            fasta.write(fasta.ungap(ref), unal)
        nseq = sum(1 for l in open(unal) if l.startswith(">"))
        for tool in a.tools.split(","):
            if (name, tool) in done:
                continue
            if tool == "magus":
                tools.TOOLS["magus"] = magus_tool(k)
            out = os.path.join(w, tool + ".fasta")
            if os.path.exists(out):
                os.remove(out)
            cwd = GCMX_PARENT if tool == "magus" else w
            res = tools.run(tool, unal, out, a.threads, w, a.timeout, os.path.join(w, tool + ".log"), cwd=cwd)
            row = {"dataset": name, "tool": tool, "nseq": nseq, "threads": a.threads, **res}
            if tool == "magus":
                row["K"] = k
            if res["status"] == "ok":
                try:
                    s = score.fastsp(true, out)
                    row.update({x: s[x] for x in ("SPFN", "SPFP", "avgErr", "TC", "LenEst", "LenRef")})
                except Exception as e:
                    row["status"] = "score failed: " + repr(e)[:200]
            print(json.dumps(row), flush=True)
            with open(a.out, "a") as f:
                f.write(json.dumps(row) + "\n")
            if not a.keep and os.path.exists(out):
                os.remove(out)
            shutil.rmtree(os.path.join(w, "magus_work"), ignore_errors=True)
            shutil.rmtree(os.path.join(w, "twilight_tmp"), ignore_errors=True)


if __name__ == "__main__":
    main()
