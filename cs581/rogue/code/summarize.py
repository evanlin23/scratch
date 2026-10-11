"""Tables for REPORT.md from results/rows_*.jsonl.  python summarize.py > ../results/SUMMARY.md"""

import glob
import json
import os
from collections import defaultdict

import numpy as np
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
rows = {}
for f in sorted(glob.glob(os.path.join(HERE, "..", "results", "rows_*.jsonl"))):
    for line in open(f):
        r = json.loads(line)
        rows[(r["instance"], r["aligner"], r["cond"])] = r  # dedupe (true/detect rows may repeat)


# detection stats are recomputed from the saved flagged sets (results/detect/INSTANCE.json)
for (i, al, c), r in list(rows.items()):
    p = os.path.join(HERE, "..", "results", "detect", i + ".json")
    if al == "detect" and os.path.exists(p):
        F = set(json.load(open(p))[c])
        R = set(json.load(open(p))["rogues"])
        r.update(nflag=len(F), TP=len(F & R), precision=len(F & R) / len(F) if F else None, recall=len(F & R) / len(R))


def setting(inst):
    return inst.rsplit("_R", 1)[0]


SETTINGS = sorted({setting(i) for i, _, _ in rows})


def paired(aligner, a, b, metric, sett, key_b=None):
    """a minus b over instances of `sett` having both; returns n, mean a, mean b, diff, W/T/L (a better), p."""
    xs, ys = [], []
    for (i, al, c), r in rows.items():
        if al != aligner or c != a or setting(i) != sett:
            continue
        rb = rows.get((i, aligner, b)) if b else r
        if rb is None:
            continue
        x = r.get(metric)
        y = rb.get(key_b or metric)
        if x is None or y is None:
            continue
        xs.append(x)
        ys.append(y)
    if not xs:
        return None
    xs, ys = np.array(xs) * 100, np.array(ys) * 100
    d = xs - ys
    w, t, l = int((d < -1e-9).sum()), int((abs(d) <= 1e-9).sum()), int((d > 1e-9).sum())
    try:
        p = wilcoxon(xs, ys).pvalue if (abs(d) > 1e-9).any() else 1.0
    except ValueError:
        p = float("nan")
    return len(xs), xs.mean(), ys.mean(), d.mean(), "%d/%d/%d" % (w, t, l), p


def line(label, res):
    if res is None:
        return None
    n, a, b, d, wtl, p = res
    return "| %s | %d | %.2f | %.2f | %+.2f | %s | %.3g |" % (label, n, a, b, d, wtl, p)


HDR = "| comparison (A vs B) | n | A | B | A-B | W/T/L (A better) | Wilcoxon p |\n|---|---|---|---|---|---|---|"
out = []
for sett in SETTINGS:
    out.append("\n## %s\n" % sett)
    t = paired("true", "true", "true-C", "FN", sett)
    if t:
        out.append("True alignment, tree FN%% on non-rogues: with rogues %.2f, without %.2f (n=%d, W/T/L %s, p=%.3g)\n"
                   % (t[1], t[2], t[0], t[4], t[5]))
    for d in ("ts", "pd", "hmm"):
        rs = [r for (i, al, c), r in rows.items() if al == "detect" and c == d and setting(i) == sett]
        if rs:
            out.append("Detector %-3s: n=%d, flagged %.1f, precision %.2f, recall %.2f, wall %.0fs" % (
                d, len(rs), np.mean([r["nflag"] for r in rs]),
                np.mean([r["precision"] if r["precision"] is not None else 0 for r in rs]),
                np.mean([r["recall"] for r in rs]), np.mean([r["wall"] for r in rs])))
    for X in ("mafft", "magus", "pasta"):
        if not any(al == X and setting(i) == sett for i, al, _ in rows):
            continue
        out.append("\n### %s — alignment error (avg of SPFN, SPFP; %%) on non-rogue taxa\n" % X)
        out.append(HDR)
        for a, b, lab in (("oracle", "all", "oracle-removed vs all"),
                          ("ts", "all", "TreeShrink filter vs all"), ("pd", "all", "p-dist filter vs all"),
                          ("hmm", "all", "HMM filter vs all"), ("pd", "oracle", "p-dist filter vs oracle"),
                          ("hmm", "oracle", "HMM filter vs oracle"), ("ts", "oracle", "TreeShrink vs oracle")):
            l = line(lab, paired(X, a, b, "avgErr", sett))
            if l:
                out.append(l)
        l = line("CONTROL: 50 random non-rogues dropped vs all (on remaining)", paired(X, "rand-drop", None, "avgErr", sett, key_b="base_avgErr"))
        if l:
            out.append(l)
        for m in ("SPFN", "SPFP"):
            l = line("oracle vs all [%s]" % m, paired(X, "oracle", "all", m, sett))
            if l:
                out.append(l)
        out.append("\n### %s — tree FN (%%) on non-rogue taxa (FastTree)\n" % X)
        out.append(HDR)
        for a, b, lab in (("oracle", "all", "oracle-removed vs all"),
                          ("all-treeC", "all", "rogue rows dropped before FastTree vs all"),
                          ("oracle", "all-treeC", "oracle vs rogue rows dropped (alignment effect)"),
                          ("oracle+add", "all", "oracle + HMM add-back vs all"),
                          ("oracle+add", "oracle", "oracle + add-back vs oracle"),
                          ("ts", "all", "TreeShrink filter + add-back vs all"),
                          ("pd", "all", "p-dist filter + add-back vs all"),
                          ("hmm", "all", "HMM filter + add-back vs all")):
            l = line(lab, paired(X, a, b, "FN", sett))
            if l:
                out.append(l)
        l = line("CONTROL: 50 random non-rogues dropped vs all (on remaining)", paired(X, "rand-drop", None, "FN", sett, key_b="base_FN"))
        if l:
            out.append(l)
        for d in ("ts", "pd", "hmm"):
            l = line("%s-drop vs all (on C minus flagged)" % d, paired(X, d + "-drop", None, "FN", sett, key_b="base_FN"))
            if l:
                out.append(l)
        walls = defaultdict(list)
        for (i, al, c), r in rows.items():
            if al == X and setting(i) == sett and "wall" in r:
                walls[c].append(r["wall"])
        out.append("\nwall-clock (s, alignment step incl. add-back): " + ", ".join(
            "%s %.0f" % (c, np.mean(v)) for c, v in sorted(walls.items()) if not c.endswith("drop")))
print("# Rogue pilot: summary tables (generated by code/summarize.py)\n")
print("A-B < 0 means A has lower error. W/T/L counts replicates where A is better/tied/worse.")
print("\n".join(out))
