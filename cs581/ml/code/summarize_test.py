"""Held-out test summary: per dataset and alignment, mean FN, CPU and log-likelihood, with paired
comparisons (mean difference, wins/ties/losses, two-sided Wilcoxon signed-rank p) of every method
against RAxML-NG (one parsimony start) on the same replicates.

    python summarize_test.py [results/test.jsonl ...]  > results/test_tables.md

For Park1000M1HF the published Park et al. 2021 trees (IDB-7008049) are added as extra rows:
"published RAxML-NG (20 starts, 24 h cap)" and "published IQ-TREE 2 (default)".
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import treeerr  # noqa: E402

R = os.path.join(HERE, "..", "results")
PARK = "/opt/data/park2021/1000M1_HF"
PUB = {"published RAxML-NG (20 starts, 24 h cap)": "RAxML-ng/R{r}/raxmlng-result.raxml.lastTree.TMP",
       "published IQ-TREE 2 (default)": "IQTree/R{r}/iqtree-result.treefile",
       "published GTM (IQ-TREE start)": "GTM/500/iqtree/R{r}/branch_length."}
ORDER = ["fasttree", "iqtree_fast", "raxmlng_fastmode", "raxmlng", "frag_constr_fasttree_c0.8",
         "frag_fpolish_fasttree_c0.8"] + list(PUB)
LABEL = {"fasttree": "FastTree 2", "iqtree_fast": "IQ-TREE 3 --fast", "raxmlng_fastmode": "RAxML-NG --fast",
         "raxmlng": "RAxML-NG (1 start)", "frag_constr_fasttree_c0.8": "**ours-fast** (constrained)",
         "frag_fpolish_fasttree_c0.8": "**ours-accurate** (+ fast polish)"}
DESC = {"Park1000M1HF": "1000M1-HF, Park et al. 2021 published inputs (simulated)",
        "1000M2HF": "1000M2-HF (simulated, our fragmentation)",
        "RNASimHF": "RNASim1K-HF (simulated, our fragmentation)",
        "16SMHF": "16S.M, real 16S rRNA, random fragments; reference = CRW tree (47% resolved)",
        "16SMHFamp": "16S.M, real 16S rRNA, amplicon-like fragments; reference = CRW tree"}


def main():
    files = sys.argv[1:] or [os.path.join(R, "test.jsonl")]
    rows = [json.loads(l) for f in files for l in open(f)]
    # Park1000M1HF baselines run earlier (FastTree, IQ-TREE --fast, RAxML-NG 1 start), same settings
    for f in ("baseline.jsonl", "hf.jsonl"):
        rows += [r for r in map(json.loads, open(os.path.join(R, f))) if r["dataset"] == "Park1000M1HF"
                 and r["method"] in ("fasttree", "iqtree_fast", "raxmlng")]
    cell = defaultdict(dict)  # (ds, aln) -> method -> rep -> row
    for r in rows:
        cell[(r["dataset"], r["aln"])].setdefault(r["method"], {})[str(r["rep"])] = r
    for rep in range(5):  # published Park trees, true alignment only
        true = os.path.join("/opt/data/Datasets/ROSE/1000M1", "R%d" % rep, "rose.tt")
        for m, pat in PUB.items():
            p = os.path.join(PARK, pat.format(r=rep))
            if os.path.exists(p) and ("Park1000M1HF", "true_align") in cell:
                cell[("Park1000M1HF", "true_align")].setdefault(m, {})[str(rep)] = {
                    "fn_rate": treeerr.error(true, p)["fn_rate"]}
    for (ds, aln), meth in sorted(cell.items()):
        ref = meth.get("raxmlng", {})
        print("\n### %s — %s alignment\n" % (DESC.get(ds, ds), "true" if aln == "true_align" else aln.upper()))
        print("| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |")
        print("|---|---|---|---|---|---|---|")
        for m in ORDER:
            d = meth.get(m)
            if not d:
                continue
            fn = [v["fn_rate"] for v in d.values()]
            cpu = [v.get("cpu_seconds") or v.get("seconds") for v in d.values() if (v.get("cpu_seconds") or v.get("seconds"))]
            keys = [k for k in d if k in ref]
            diff = [100 * (d[k]["fn_rate"] - ref[k]["fn_rate"]) for k in keys] if m != "raxmlng" else []
            wtl = (sum(x < -1e-9 for x in diff), sum(abs(x) <= 1e-9 for x in diff), sum(x > 1e-9 for x in diff))
            p = "-"
            if len(diff) >= 3 and any(abs(x) > 1e-9 for x in diff):
                p = "%.2f" % wilcoxon(diff).pvalue
            dl = [d[k]["lnl_tool"] - ref[k]["lnl_tool"] for k in keys
                  if d[k].get("lnl_tool") is not None and ref[k].get("lnl_tool") is not None
                  and m.startswith(("raxmlng", "frag"))]
            print("| %s | %d | %.1f%% | %s | %s | %s | %s |" % (
                LABEL.get(m, m), len(fn), 100 * np.mean(fn),
                "%+.1f (%d/%d/%d)" % (np.mean(diff), *wtl) if diff else "–", p,
                "%.0f" % (np.mean(cpu) / 60) if cpu else "–", "%+.1f" % np.mean(dl) if dl else "–"))
        if aln == "upp":
            ts = [json.load(open(os.path.join("/opt/data/mlcache", ds, "R%s" % k, "upp_time.json")))["cpu_seconds"]
                  for k in ref if os.path.exists(os.path.join("/opt/data/mlcache", ds, "R%s" % k, "upp_time.json"))]
            if ts:
                print("\nUPP alignment cost (not included above): %.0f CPU min per replicate." % (np.mean(ts) / 60))
    pooled(cell, ["fasttree", "iqtree_fast", "raxmlng_fastmode", "frag_constr_fasttree_c0.8",
                  "frag_fpolish_fasttree_c0.8"])


def pooled(rows_by_cell, methods):
    """Pooled paired comparison over every held-out replicate (all datasets and alignments)."""
    print("\n### Pooled over all held-out replicates (paired vs RAxML-NG, 1 start)\n")
    print("| method | n pairs | mean Δ FN | W/T/L | Wilcoxon p | mean CPU ratio vs RAxML-NG |")
    print("|---|---|---|---|---|---|")
    for m in methods:
        diff, ratio = [], []
        for meth in rows_by_cell.values():
            ref, d = meth.get("raxmlng", {}), meth.get(m, {})
            for k in d:
                if k in ref:
                    diff.append(100 * (d[k]["fn_rate"] - ref[k]["fn_rate"]))
                    if d[k].get("cpu_seconds") and ref[k].get("cpu_seconds"):
                        ratio.append(d[k]["cpu_seconds"] / ref[k]["cpu_seconds"])
        if len(diff) < 3:
            continue
        wtl = (sum(x < -1e-9 for x in diff), sum(abs(x) <= 1e-9 for x in diff), sum(x > 1e-9 for x in diff))
        print("| %s | %d | %+.2f | %d/%d/%d | %.3g | %.2f |" % (LABEL.get(m, m), len(diff), np.mean(diff), *wtl,
              wilcoxon(diff).pvalue, np.mean(ratio) if ratio else float("nan")))


if __name__ == "__main__":
    main()
