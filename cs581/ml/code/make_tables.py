"""Regenerate every results table used in REPORT.md:  python make_tables.py > ../results/tables.md"""
import json
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "results")


def load(*names):
    rows = []
    for n in names:
        p = os.path.join(R, n)
        if os.path.exists(p):
            rows += [json.loads(l) for l in open(p)]
    return rows


def sh(*args):
    return subprocess.run([sys.executable] + list(args), capture_output=True, text=True, cwd=HERE).stdout


def paired(rows, ds, ref, methods, title):
    """Per-replicate paired comparison against method `ref` (same dataset/rep/alignment)."""
    by = defaultdict(dict)
    for r in rows:
        if r["dataset"] == ds:
            by[(str(r["rep"]), r["aln"])][r["method"]] = r
    print("\n### %s\n" % title)
    print("| method | n | mean FN | FN vs %s (mean paired diff, wins/ties/losses) | mean CPU min | mean lnL - %s lnL |" % (ref, ref))
    print("|---|---|---|---|---|---|")
    for m in methods:
        fn, dfn, cpu, dl, w = [], [], [], [], [0, 0, 0]
        for key, d in sorted(by.items()):
            if m not in d:
                continue
            fn.append(d[m]["fn_rate"])
            c = d[m].get("cpu_seconds") or d[m].get("seconds")
            if m.startswith("raxmlng_from_") or m == "raxmlng_ft":  # add the start tree's cost
                src = "fasttree" if m == "raxmlng_ft" else m.split("_from_")[1]
                if src in d:
                    c += d[src].get("cpu_seconds") or d[src]["seconds"]
            cpu.append(c / 60)
            if ref in d and m != ref:
                x = d[m]["fn_rate"] - d[ref]["fn_rate"]
                dfn.append(x)
                w[0 if x < -1e-9 else (1 if abs(x) <= 1e-9 else 2)] += 1
                if d[m].get("lnl_tool") is not None and d[ref].get("lnl_tool") is not None and m.startswith(("raxmlng", "frag")):
                    dl.append(d[m]["lnl_tool"] - d[ref]["lnl_tool"])
        if not fn:
            continue
        print("| %s | %d | %.1f%% | %s | %.1f | %s |" % (
            m, len(fn), 100 * np.mean(fn),
            "%+.1f pts (%d/%d/%d)" % (100 * np.mean(dfn), *w) if dfn else "-",
            np.mean(cpu), "%+.1f" % np.mean(dl) if dl else "-"))


print("## Tree error (FN) on the MAGUS-paper data, mean over R0-R2\n")
print(sh("summarize.py", os.path.join(R, "baseline.jsonl"), os.path.join(R, "masking.jsonl")))
print("\n## CPU seconds\n")
print(sh("summarize.py", os.path.join(R, "baseline.jsonl"), "--metric", "cpu_seconds"))
print("\n## Validation: published Park et al. 2021 trees rescored\n")
print(open(os.path.join(R, "validate_park.md")).read())
print("\n## Validation: PASTA's published tree vs our FastTree on the PASTA alignment\n")
print(sh("compare_pasta_tre.py", os.path.join(R, "baseline.jsonl")))
hf = load("baseline.jsonl", "hf.jsonl", "frag.jsonl")
paired(hf, "Park1000M1HF", "raxmlng",
       ["fasttree", "iqtree_fast", "raxmlng", "raxmlng_from_iqtree_fast", "frag_constr_fasttree", "frag_polish_fasttree"],
       "Pilot: 1000M1-HF (true alignment, half the sequences fragmentary), paired vs RAxML-NG (1 parsimony start)")
st = load("baseline.jsonl", "starttree.jsonl", "timing_probe.jsonl")
for ds in ("1000M2", "1000M3", "1000L1", "RNASim"):
    paired([r for r in st if r["aln"] == "gcm" and str(r["rep"]) == "0"], ds, "raxmlng",
           ["fasttree", "iqtree_fast", "raxmlng", "raxmlng_ft"], "Pilot: start trees, %s R0, MAGUS alignment" % ds)
mk = [r for r in load("baseline.jsonl", "masking.jsonl") if r["method"] == "fasttree"]
fn = {(r["dataset"], str(r["rep"]), r["aln"]): r["fn_rate"] for r in mk}
print("\n### Pilot: column masking / alignment ensembles (FastTree), FN change vs the unmasked MAGUS alignment\n")
variants = ["gcm_mask50", "gcm_mask70", "gcm_oracle50", "gcm_oracle70", "gcm_pasta", "pasta_align", "true_align"]
dss = ["1000M2", "1000M3", "1000L1", "RNASim"]
print("| alignment given to FastTree | " + " | ".join(dss) + " | all 12 (wins/ties/losses) |")
print("|---|" + "---|" * (len(dss) + 1))
for v in variants:
    cells, alld = [], []
    for ds in dss:
        d = [fn[(ds, rep, v)] - fn[(ds, rep, "gcm")] for rep in "012" if (ds, rep, v) in fn and (ds, rep, "gcm") in fn]
        alld += d
        cells.append("%+.1f" % (100 * np.mean(d)) if d else "-")
    w = (sum(x < -1e-9 for x in alld), sum(abs(x) <= 1e-9 for x in alld), sum(x > 1e-9 for x in alld))
    print("| %s | %s | %+.1f pts (%d/%d/%d) |" % (v, " | ".join(cells), 100 * np.mean(alld), *w))
