# AI-assisted (Claude), exploration code for CS581 project
"""Pool FastTree nRF rows: this branch (/opt/work/gcmvote/reps/*/trees.jsonl) + helper branches
claude/cs581-gcmvote-* (cs581/gcmvote/results_*/trees.jsonl). Writes ../results/trees.md and
../results/helpers/<branch>.trees.jsonl.

Method names are normalised: vote.py variants -> vote_<variant>; the masked alignments -> vote_hard+mask,
magus+mask; the helpers' gcmtrees baselines keep their names (recipe = wsoft0.03:linsi&fftns2#es4,
hard = linsi&fftns2-op3 -> renamed fftns2-hardfilter).
"""
import glob
import json
import os
import subprocess

import numpy as np
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
MINE = {"magus": "magus", "es4": "es4", "true": "true", "hard": "vote_hard", "hard-bb": "vote_hard-bb",
        "hard_mask_masked": "vote_hard+mask", "magus_mask_masked": "magus+mask"}


def helper_rows():
    git = lambda *a: subprocess.run(["git", "-C", REPO] + list(a), capture_output=True, text=True)
    refs = git("ls-remote", "origin", "refs/heads/claude/cs581-gcmvote-*").stdout.split()
    rows = []
    os.makedirs(os.path.join(OUT, "helpers"), exist_ok=True)
    for ref in refs[1::2]:
        br = ref.replace("refs/heads/", "")
        git("fetch", "-q", "origin", "{0}:refs/remotes/origin/{0}".format(br))
        files = [f for f in git("ls-tree", "-r", "--name-only", "origin/" + br).stdout.split()
                 if f.startswith("cs581/gcmvote/results_") and f.endswith("trees.jsonl")]
        for f in files:
            txt = git("show", "origin/{}:{}".format(br, f)).stdout
            open(os.path.join(OUT, "helpers", br.split("/")[-1] + "." + os.path.basename(f)), "w").write(txt)
            for l in txt.splitlines():
                if l.strip():
                    d = json.loads(l)
                    m = d["method"]
                    m = {"hard": "fftns2-hardfilter", "vote_hard_mask": "vote_hard+mask"}.get(m, m)
                    rows.append((d["dataset"], m, 100 * d["RF"], br.split("-")[-1]))
    return rows


def main():
    T, src = {}, {}
    for f in glob.glob("/opt/work/gcmvote/reps/*/trees.jsonl"):
        for l in open(f):
            d = json.loads(l)
            T[(d["dataset"], MINE.get(d["method"], d["method"]))] = 100 * d["RF"]
            src[d["dataset"]] = "this"
    for ds, m, rf, br in helper_rows():
        if src.get(ds, br) != br:  # dataset already from another source
            continue
        src[ds] = br
        T[(ds, m)] = rf
    ds = sorted({d for d, _ in T}, key=lambda x: (x.split("_R")[0], int(x.split("_R")[1])))
    ms = ["true", "magus", "es4", "vote_hard-bb", "vote_hard", "vote_hard+mask", "magus+mask", "vote_soft",
          "vote_soft4", "recipe", "es3"]
    ms = [m for m in ms if any((d, m) in T for d in ds)]
    out = ["# Trees: FastTree -lg -gamma nRF (%) vs the true tree, pooled\n",
           "Rows from this branch (`this`) and from helper branches claude/cs581-gcmvote-<x> (one MAGUS draw each).\n",
           "| dataset | source | " + " | ".join(ms) + " |", "|---|---|" + "---|" * len(ms)]
    for d in ds:
        out.append("| {} | {} | {} |".format(d, src[d], " | ".join("{:.2f}".format(T[(d, m)]) if (d, m) in T else "–" for m in ms)))
    out += ["", "| method | Δ nRF vs magus: mean (W/T/L at 0.1, two-sided Wilcoxon p), n | SIMHIGH R1–R4 only |", "|---|---|---|"]
    for m in ms:
        if m == "magus":
            continue
        cells = []
        for sel in (lambda d: True, lambda d: d.startswith("SIMHIGH_R") and int(d.split("_R")[1]) <= 4):
            dd = [T[(d, m)] - T[(d, "magus")] for d in ds if sel(d) and (d, m) in T and (d, "magus") in T]
            if not dd:
                cells.append("–")
                continue
            w = sum(x < -0.1 for x in dd); l = sum(x > 0.1 for x in dd)
            p = wilcoxon(dd).pvalue if len(dd) > 1 and any(dd) else float("nan")
            cells.append("{:+.2f} ({}/{}/{}, p = {:.3f}), n = {}".format(np.mean(dd), w, len(dd) - w - l, l, p, len(dd)))
        out.append("| {} | {} | {} |".format(m, *cells))
    open(os.path.join(OUT, "trees.md"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
