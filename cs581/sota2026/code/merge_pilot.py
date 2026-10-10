"""Pilot of "MAGUS with a different base method": GCM merge with the base aligner swapped.

    python3 merge_pilot.py OUT.jsonl REP[,REP...] TOOL[,TOOL...] [--backbones cached|tool] [--threads 4]

For each replicate with a cached MAGUS run (cs581/experiments/runs/<rep>/inputs.tar.xz):
  1. realign each of MAGUS's 25 decomposition subsets with TOOL (same sequences, from scratch);
  2. backbones: "cached" = MAGUS's own MAFFT L-INS-i backbones (only the subset aligner changes),
     "tool" = the same 10 backbone sequence sets realigned with TOOL (TOOL is the whole base method);
  3. run MAGUS's GCM merge only (gcmx.run_magus -s ... -b ..., paper's merge flags);
  4. score with FastSP against the full reference.
TOOL "cached-linsi" reuses MAGUS's own subset alignments: it is the control and reproduces MAGUS.
Times: subset+backbone alignment wall-clock (sequential, `threads` threads each) and merge wall-clock.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import subsets as S  # noqa: E402  (reference(), cached(), fasta, score, tools)

fasta, score, tools = S.fasta, S.score, S.tools
CODE = os.path.join(S.REPO, "cs581", "code")
MERGE_FLAGS = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
               "-f", "4"]


def realign(tool, src_dir, dst_dir, threads, work):
    os.makedirs(dst_dir, exist_ok=True)
    total = 0.0
    for f in sorted(os.listdir(src_dir)):
        src, dst = os.path.join(src_dir, f), os.path.join(dst_dir, f)
        if tool == "cached-linsi":
            shutil.copy(src, dst)
            continue
        unal = os.path.join(work, "u_" + f)
        fasta.write(fasta.ungap(fasta.upper(fasta.read(src))), unal)
        tw = os.path.join(work, "tw_" + f)
        os.makedirs(tw, exist_ok=True)
        r = tools.run(tool, unal, dst, threads, tw, 3600, os.path.join(work, f + ".log"))
        shutil.rmtree(tw, ignore_errors=True)
        if r["status"] != "ok":
            raise RuntimeError("{} failed on {}: {}".format(tool, f, r["status"]))
        total += r["wall"]
    return total


def main():
    p = argparse.ArgumentParser()
    p.add_argument("out")
    p.add_argument("reps")
    p.add_argument("tools")
    p.add_argument("--backbones", default="cached", choices=("cached", "tool"))
    p.add_argument("--threads", type=int, default=4)
    a = p.parse_args()
    done = set()
    if os.path.exists(a.out):
        done = {(r["rep"], r["tool"], r["backbones"]) for r in map(json.loads, open(a.out))}
    for rep in a.reps.split(","):
        inp = S.cached(rep)
        ref = fasta.upper(fasta.read(S.reference(rep)))
        for tool in a.tools.split(","):
            bbmode = "cached" if tool == "cached-linsi" else a.backbones
            if (rep, tool, bbmode) in done:
                continue
            w = os.path.join(S.WORK, "merge", rep, tool + "_" + bbmode)
            shutil.rmtree(w, ignore_errors=True)
            os.makedirs(w)
            row = {"rep": rep, "tool": tool, "backbones": bbmode}
            try:
                row["subset_wall"] = round(realign(tool, os.path.join(inp, "subalignments"), os.path.join(w, "subs"),
                                                   a.threads, w), 1)
                if bbmode == "tool":
                    row["backbone_wall"] = round(realign(tool, os.path.join(inp, "backbones"), os.path.join(w, "bbs"),
                                                         a.threads, w), 1)
                    bb = os.path.join(w, "bbs")
                else:
                    bb = os.path.join(inp, "backbones")
                # subset-level error of the realigned subsets (true alignment restricted to each subset)
                errs = []
                for f in os.listdir(os.path.join(w, "subs")):
                    est = os.path.join(w, "subs", f)
                    tr = os.path.join(w, "true_" + f)
                    fasta.write(fasta.restrict(ref, list(fasta.read(est))), tr)
                    errs.append(score.fastsp(tr, est)["avgErr"])
                row["subset_err_mean"] = sum(errs) / len(errs)
                out = os.path.join(w, "merged.fasta")
                start = time.time()
                with open(os.path.join(w, "merge.log"), "w") as log:
                    subprocess.run([sys.executable, "-m", "gcmx.run_magus", "-np", str(a.threads), "-d",
                                    os.path.join(w, "magus"), "-s", os.path.join(w, "subs"), "-b", bb, "-o", out]
                                   + MERGE_FLAGS, cwd=CODE, stdout=log, stderr=subprocess.STDOUT, check=True)
                row["merge_wall"] = round(time.time() - start, 1)
                truef = os.path.join(w, "true.fasta")
                fasta.write(ref, truef)
                s = score.fastsp(truef, out)
                row.update({k: s[k] for k in ("SPFN", "SPFP", "avgErr", "TC", "LenEst", "LenRef")})
                row["status"] = "ok"
            except Exception as e:
                row["status"] = repr(e)[:300]
            print(json.dumps(row), flush=True)
            with open(a.out, "a") as f:
                f.write(json.dumps(row) + "\n")
            shutil.rmtree(os.path.join(w, "magus"), ignore_errors=True)


if __name__ == "__main__":
    main()
