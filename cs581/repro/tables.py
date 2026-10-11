# AI-assisted (Claude), code for CS581 project
"""Core tables of the GCM evidence-filter experiment, rebuilt from the rep bank, and a check against the logged rows.

    python3 tables.py [--data cs581/repro/data] [--out cs581/repro/out]

Reads   OUT/rows/SRC/REP/VARIANT@B.json   (run_all.sh aln / bb: one rebuilt merge each, FastSP scored)
        OUT/trees/SRC/REP/METHOD.nwk      (run_all.sh trees: FastTree -lg -gamma)
        DATA/logged/...                   (fetch_bank.sh: the rows the experiment branches logged)
        DATA/bank/SRC/REP/true_tree.nwk
Writes  OUT/TABLES.md       training / held-out / extra replicates (magus, es4, vote-hard-bb), trees, backbone count
        OUT/MISMATCHES.md   every rebuilt number next to the logged number for the same MAGUS draw, plus the
                            logged summary tables recomputed from the logged rows

Error = (SPFN + SPFP) / 2 x 100 (FastSP). Delta = variant - magus (same replicate, same draw).
Summary cells: mean Delta (W/T/L, two-sided Wilcoxon signed-rank p), n; W = Delta < -band, L = Delta > +band,
band 0.05 for alignments (as gcmvote/code/summarize.py) and 0.1 for nRF (as gcmvote/code/pool_trees.py).
"""

import argparse
import glob
import json
import math
import os
import re

import numpy as np
from scipy.stats import wilcoxon

# gcmgen's fixed split (cs581/gcmgen/SPLIT.md, gcmvote PREREG.md)
TRAIN = ["BBA0101", "BBA0134", "BBA0067", "BBA0039", "SIMMOD_R1", "SIMHIGH_R1", "1000M2", "1000L1", "1000L2", "16S.M"]
HELD = ["BBA0154", "BBA0190", "SIMMOD_R2", "SIMHIGH_R2", "1000L3", "1000M3", "1000S1", "1000S2", "RNASim", "1000M4",
        "1000S3", "1000M2_R1", "1000L1_R1", "RNASim_R1"]
CORE = ["magus", "es4", "vote-hard-bb"]
LOGGED_NAME = {"vote-hard-bb": "hard-bb", "vote-hard": "hard"}  # filter_gcm.py name -> vote.py name

# Logged summary numbers that come from cs581/literature/ranking_final.md (not from a results file).
# Check-in 12 (03:45 UTC): pooled held-out on the helper draws, tie band 0.1.
CHECKIN12 = {
    "proteins": {"reps": ["SIMMOD_R2", "SIMHIGH_R2", "BBA0154", "BBA0190", "BBA0081", "BBA0117", "HF_blmb", "HF_aat",
                          "HF_Acetyltransf", "HF_PDZ"],
                 "hard-bb - magus": "-1.66 (median -1.65, 8/0/2, p = 0.084)",
                 "es4 - magus": "-1.89 (median -0.68, 6/2/2, p = 0.049)",
                 "hard-bb - es4": "+0.22 (median -0.17, 5/1/4, p = 0.92)"},
    "DNA/RNA": {"reps": ["1000L3", "1000M3", "1000M4", "1000S3", "1000M2_R1", "1000L1_R1", "1000M1_R0", "1000M1_R1",
                         "1000M1_R2"],
                "hard-bb - magus": "+0.38 (1/3/5, p = 0.074)", "es4 - magus": "-0.05",
                "hard-bb - es4": "+0.43 (p = 0.055)"},
}
# Backbone-count sensitivity (h6, h6b, h9, h9b), SP error %, B = 5 / 10 / 20: MAGUS, best fixed fraction, hard-bb
LOGGED_BB = {
    "BBA0101": ([28.79, 28.78, 28.83], [27.86, 28.62, 28.63], [27.40, 27.65, 28.29]),
    "BBA0067": ([26.40, 26.23, 26.36], [26.40, 26.11, 26.15], [26.10, 25.55, 25.62]),
    "16S.M": ([13.12, 13.03, 13.04], [13.04, 12.91, 12.95], [12.25, 12.37, 12.25]),
    "SIMHIGH_R1": ([24.57, 24.05, 24.78], [18.48, 17.88, 17.77], [18.99, 19.21, 20.23]),
    "SIMMOD_R1": ([12.88, 13.59, 13.33], [9.40, 9.54, 9.40], [10.84, 9.24, 9.38]),
    "1000M2": ([10.09, 9.00, 8.62], [10.09, 8.94, 8.57], [10.24, 9.62, 8.68]),
    "1000L1": ([8.16, 7.82, 7.52], [7.87, 7.66, 7.31], [7.72, 7.57, 7.14]),
}
# REPORT section 6 (gcmvote/results/trees.md): pooled FastTree nRF over 22 datasets, delta vs magus, band 0.1
LOGGED_TREES = {"true": "-2.81 (19/1/2, p = 0.000), n = 22", "es4": "-0.13 (12/1/9, p = 0.355), n = 22",
                "vote-hard-bb": "-0.23 (10/3/9, p = 0.614), n = 22"}


# ---------------------------------------------------------------------------------------------- names and groups

def dataset_of(rep):
    """Bank replicate name -> dataset name used in the tables ('1000L3_R0' -> '1000L3', 'BBA0101_b20' -> 'BBA0101',
    '16S.M_R0_B10' -> '16S.M'); replicate numbers other than R0 are kept (1000M1_R2, SIMHIGH_R7)."""
    d = re.sub(r"_(B(5|10|20)|b20)$", "", rep)
    if d.endswith("_R0") and (d[:-3] in TRAIN + HELD or d.startswith("16S.")):
        d = d[:-3]
    return d


def group_of(ds):
    if ds.startswith(("BBA", "SIM", "HF_")):
        return "proteins"
    if ds.startswith("16S."):
        return "real rRNA"
    return "simulated DNA/RNA"


GROUPS = ["proteins", "simulated DNA/RNA", "real rRNA"]


# ---------------------------------------------------------------------------------------------- loading

def load_rebuilt(out):
    rows = {}
    for f in glob.glob(os.path.join(out, "rows", "*", "*", "*.json")):
        r = json.load(open(f))
        rows[(r["src"], r["rep"], r["variant"], r["B"])] = r
    return rows


def load_logged(data):
    """Every logged alignment row, normalised to dict(src, rep, impl, variant, B, avgErr, SPFN, SPFP, file).
    impl = 'vote' for gcmvote/code/vote.py rows (the code filter_gcm.py copies), 'gg' for gcmgen/code/gg.py rows
    (MAGUS's merge = 'linsi', es filter = 'linsi#esK'; gg runs MAGUS with its own dict graph, so its es rows can
    differ slightly from vote.py's)."""
    rows = []
    for f in sorted(glob.glob(os.path.join(data, "logged", "**", "*.jsonl"), recursive=True)):
        src = os.path.relpath(f, os.path.join(data, "logged")).split(os.sep)[0]
        for line in open(f):
            if not line.strip():
                continue
            r = json.loads(line)
            if "avgErr" not in r or "variant" not in r:
                continue
            rep = r.get("rep") or r.get("dataset")
            v = r["variant"]
            B = int(r.get("B") or r.get("nbb") or 10)
            if r.get("gg_variant"):
                v = r["gg_variant"]
            if v.startswith("vote:"):
                impl, v = "vote", v[5:]
            elif v == "linsi":
                impl, v = "gg", "magus"
            elif re.fullmatch(r"linsi#es\d+", v):
                impl, v = "gg", v.split("#")[1]
            elif re.fullmatch(r"magus|es\d+|frac[0-9.]+|hard(-bb)?|hard[0-9.]+|soft\d*(-bb)?|hard\+mask", v):
                impl = "vote"
            else:
                continue
            if r.get("src") == "gg":
                impl = "gg"
            rows.append(dict(src=src, rep=rep, impl=impl, variant=v, B=B, avgErr=r["avgErr"], SPFN=r.get("SPFN"),
                             SPFP=r.get("SPFP"), file=os.path.relpath(f, data)))
    return rows


def logged_match(logged, src, rep, variant, B):
    """The logged row(s) for the same draw (same source branch) and the same rule. Returns (row, how)."""
    names = {rep, re.sub(r"_(B(5|10|20)|b20)$", "", rep)}
    v = LOGGED_NAME.get(variant, variant)
    cand = [r for r in logged if r["src"] == src and r["rep"] in names and r["B"] == B]
    exact = [r for r in cand if r["impl"] == "vote" and r["variant"] == v]
    if exact:
        return exact[0], "vote.py"
    if v.startswith("frac"):  # fixed fraction F at B = es with K = ceil(F * B)
        v = "es{}".format(math.ceil(float(v[4:]) * B - 1e-9))
        exact = [r for r in cand if r["impl"] == "vote" and r["variant"] == v]
        if exact:
            return exact[0], "vote.py"
    gg = [r for r in cand if r["impl"] == "gg" and r["variant"] == v]
    if gg:
        return gg[0], "gg.py"
    return None, None


# ---------------------------------------------------------------------------------------------- statistics

def summary(ds, band=0.05):
    ds = [d for d in ds if d is not None]
    if not ds:
        return "–"
    w = sum(d < -band for d in ds)
    l = sum(d > band for d in ds)
    t = len(ds) - w - l
    p = wilcoxon(ds).pvalue if len(ds) > 1 and any(abs(d) > 0 for d in ds) else float("nan")
    return "{:+.2f} ({}/{}/{}, p = {:.3f}), n = {}".format(np.mean(ds), w, t, l, p, len(ds))


def table(lines, title, reps, err, variants, base="magus", band=0.05, groups=None):
    """Per-replicate Delta table plus one summary row per group. err(rep, variant) -> error % or None."""
    lines += ["### " + title, "", "| replicate | source | {} % | {} |".format(base, " | ".join(
        "{} − {}".format(v, base) for v in variants)), "|---|---|---|" + "---|" * len(variants)]
    kept = []
    for ds, src in reps:
        b = err(ds, base)
        if b is None:
            continue
        kept.append(ds)
        cells = ["{:+.2f}".format(err(ds, v) - b) if err(ds, v) is not None else "–" for v in variants]
        lines.append("| {} | {} | {:.2f} | {} |".format(ds, src, b, " | ".join(cells)))
    groups = groups or [(g, lambda d, g=g: group_of(d) == g) for g in GROUPS] + \
        [("DNA/RNA (simulated + rRNA)", lambda d: group_of(d) != "proteins"), ("all", lambda d: True)]
    for lab, sel in groups:
        ds = [d for d in kept if sel(d)]
        if ds:
            cells = [summary([err(d, v) - err(d, base) if err(d, v) is not None else None for d in ds], band)
                     for v in variants]
            lines.append("| **{}** | | | {} |".format(lab, " | ".join(cells)))
    lines.append("")


# ---------------------------------------------------------------------------------------------- trees

def splits(path):
    """Set of non-trivial bipartitions of an unrooted newick tree, each as a frozenset of leaf names on the side
    without the alphabetically first leaf, plus the leaf set. Branch lengths / support labels are ignored."""
    text = open(path).read().strip().rstrip(";")
    tokens = re.findall(r"\(|\)|,|:[^,()]*|[^,():]+", text)
    stack, clades, after_close = [[]], [], False
    for tok in tokens:
        if tok == "(":
            stack.append([])
            after_close = False
        elif tok == ")":
            members = stack.pop()
            clade = frozenset().union(*members)
            clades.append(clade)
            stack[-1].append(clade)
            after_close = True
        elif tok == ",":
            after_close = False
        elif tok.startswith(":"):
            continue
        elif not after_close:  # a leaf name (a label right after ')' is a support value)
            stack[-1].append(frozenset([tok.strip()]))
    leaves = frozenset().union(*stack[0])
    first = min(leaves)
    out = set()
    for c in clades:
        side = leaves - c if first in c else c
        if 1 < len(side) < len(leaves) - 1:
            out.add(frozenset(side))
    return out, leaves


def nrf(true_tree, est_tree):
    """(FN + FP) / (2 (n - 3)), as cs581/protbench/code/trees.py computes it with DendroPy."""
    t, lt = splits(true_tree)
    e, le = splits(est_tree)
    assert lt == le, "trees have different leaves"
    return round((len(t - e) + len(e - t)) / (2 * (len(lt) - 3)), 4)


# ---------------------------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--data", default=os.path.join(here, "data"))
    ap.add_argument("--out", default=os.path.join(here, "out"))
    a = ap.parse_args()
    rebuilt = load_rebuilt(a.out)
    logged = load_logged(a.data)
    T, M = [], []  # TABLES.md, MISMATCHES.md lines

    # which bank replicate stands for which dataset: gcmvote helper draws first (their vote rows were logged),
    # gcmtrees draws only for datasets no gcmvote branch banked
    reps_at_b10 = sorted({(s, r) for (s, r, v, B) in rebuilt if B == 10})
    pick = {}
    for src, rep in sorted(reps_at_b10, key=lambda x: (not x[0].startswith("gcmvote"), x)):
        if re.search(r"_B(5|20)$", rep):
            continue
        pick.setdefault(dataset_of(rep), (src, rep))

    def err(ds, v, B=10):
        src, rep = pick.get(ds, (None, None))
        r = rebuilt.get((src, rep, v, B))
        return None if r is None else 100 * r["avgErr"]

    def listed(names):
        return [(d, pick[d][0]) for d in names if d in pick]

    T += ["<!-- AI-assisted (Claude), code for CS581 project; generated by tables.py -->", "",
          "# Rebuilt tables (from the rep bank)", "",
          "Error = (SPFN + SPFP) / 2 × 100. Δ = variant − magus on the same MAGUS draw. Summary: mean Δ "
          "(W/T/L at ±0.05, two-sided Wilcoxon p), n. `source` = the branch whose banked MAGUS draw was used.", ""]
    missing_train = [d for d in TRAIN if d not in pick]
    table(T, "Training split (B = 10)", listed(TRAIN), err, ["es4", "vote-hard-bb"])
    T += ["Training datasets with no banked draw: {}.".format(", ".join(missing_train) or "none"), ""]
    table(T, "Held-out split, pre-registered 14 (B = 10)", listed(HELD), err, ["es4", "vote-hard-bb"])
    table(T, "Held-out: vote-hard-bb vs es4", listed(HELD), err, ["vote-hard-bb"], base="es4")
    extra = sorted(d for d in pick if d not in TRAIN + HELD)
    table(T, "All other banked replicates (B = 10)", listed(extra), err, ["es4", "vote-hard-bb"])
    every = sorted(pick)
    table(T, "Every banked dataset, by data type (B = 10)", listed(every), err, ["es4", "vote-hard-bb"],
          groups=[(g, lambda d, g=g: group_of(d) == g) for g in GROUPS] + [("all", lambda d: True)])

    # check-in 12 pooled sets (tie band 0.1)
    T += ["### Logged pooled held-out of check-in 12 (ranking_final.md), recomputed on the same draws", "",
          "| set | comparison | rebuilt (band 0.1) | logged |", "|---|---|---|---|"]
    for lab, spec in CHECKIN12.items():
        ds = [d for d in spec["reps"] if d in pick]
        for comp, (v, b) in {"hard-bb - magus": ("vote-hard-bb", "magus"), "es4 - magus": ("es4", "magus"),
                             "hard-bb - es4": ("vote-hard-bb", "es4")}.items():
            dd = [err(d, v) - err(d, b) for d in ds if err(d, v) is not None and err(d, b) is not None]
            s = summary(dd, 0.1)
            if dd:
                s = s.replace(" (", " (median {:+.2f}, ".format(float(np.median(dd))), 1)
            T.append("| {} (n = {}) | {} | {} | {} |".format(lab, len(ds), comp, s, spec[comp]))
    T.append("")

    # backbone count
    T += ["### Backbone count B (same subsets; B = 5 is backbones 1–5, B = 20 adds 10 new L-INS-i backbones)", "",
          "SP error %, B = 5 / 10 / 20. Best fixed fraction = min over frac0.2–0.5 (k ≥ ⌈F·B⌉) in hindsight. "
          "Logged = ranking_final.md check-in 12.", "",
          "| dataset | source | MAGUS | best fixed fraction | vote-hard-bb | logged MAGUS | logged best fixed | logged hard-bb |",
          "|---|---|---|---|---|---|---|---|"]
    bb_reps = {}
    for (src, rep, v, B) in rebuilt:
        if src in ("gcmvote-h6", "gcmvote-h6b", "gcmvote-h9", "gcmvote-h9b"):
            bb_reps.setdefault(dataset_of(rep), {}).setdefault(B, (src, rep))
    for ds in sorted(bb_reps):
        def e(v, B):
            src, rep = bb_reps[ds].get(B, (None, None))
            r = rebuilt.get((src, rep, v, B))
            return None if r is None else 100 * r["avgErr"]
        fmt = lambda xs: " / ".join("–" if x is None else "{:.2f}".format(x) for x in xs)
        mag = [e("magus", B) for B in (5, 10, 20)]
        fixed = [min([x for x in (e("frac{}".format(F), B) for F in (0.2, 0.3, 0.4, 0.5)) if x is not None], default=None)
                 for B in (5, 10, 20)]
        hbb = [e("vote-hard-bb", B) for B in (5, 10, 20)]
        lg = LOGGED_BB.get(ds, ([None] * 3,) * 3)
        T.append("| {} | {} | {} | {} | {} | {} | {} | {} |".format(ds, bb_reps[ds][10][0] if 10 in bb_reps[ds] else "",
                 fmt(mag), fmt(fixed), fmt(hbb), fmt(lg[0]), fmt(lg[1]), fmt(lg[2])))
    T.append("")

    # trees
    tree_rows = {}
    for f in glob.glob(os.path.join(a.out, "trees", "*", "*", "*.nwk")):
        src, rep, m = f.split(os.sep)[-3], f.split(os.sep)[-2], os.path.basename(f)[:-4]
        tt = os.path.join(a.data, "bank", src, rep, "true_tree.nwk")
        tree_rows[(src, rep, m)] = 100 * nrf(tt, f)
    trees_by_ds = {}
    for (src, rep, m), v in tree_rows.items():
        trees_by_ds.setdefault(rep, {"src": src})[m] = v
    if trees_by_ds:
        ms = ["true", "magus", "es4", "vote-hard-bb"]
        T += ["### Trees: FastTree -lg -gamma nRF % vs the true tree (B = 10 alignments)", "",
              "| dataset | source | " + " | ".join(ms) + " |", "|---|---|" + "---|" * len(ms)]
        order = sorted(trees_by_ds, key=lambda x: (x.split("_R")[0], int(x.split("_R")[1])))
        for rep in order:
            d = trees_by_ds[rep]
            T.append("| {} | {} | {} |".format(rep, d["src"], " | ".join(
                "{:.2f}".format(d[m]) if m in d else "–" for m in ms)))
        T += ["", "| method | Δ nRF vs magus, all rebuilt (band 0.1) | same, SIMHIGH only | logged (22 datasets) |",
              "|---|---|---|---|"]
        for m in ["true", "es4", "vote-hard-bb"]:
            def dl(sel):
                return [d[m] - d["magus"] for r, d in trees_by_ds.items() if sel(r) and m in d and "magus" in d]
            T.append("| {} | {} | {} | {} |".format(m, summary(dl(lambda r: True), 0.1),
                                                    summary(dl(lambda r: r.startswith("SIMHIGH")), 0.1),
                                                    LOGGED_TREES.get(m, "")))
        T.append("")

    # ------------------------------------------------------------------ MISMATCHES.md
    M += ["<!-- AI-assisted (Claude), code for CS581 project; generated by tables.py -->", "",
          "# Rebuilt vs logged", "",
          "## 1. Every rebuilt merge vs the logged row of the same MAGUS draw", "",
          "`exact` = |Δ avgErr| < 1e-12. `vote.py` rows were logged by the code filter_gcm.py copies; `gg.py` rows "
          "by gcmgen's gg.py, whose MAGUS run writes the graph file in a different line order, so small "
          "differences are expected for its es/fixed-fraction rows.", ""]
    counts = {}
    detail = []
    for key in sorted(rebuilt):
        src, rep, v, B = key
        r = rebuilt[key]
        lg, how = logged_match(logged, src, rep, v, B)
        if lg is None:
            status = "no logged row"
        else:
            diff = 100 * (r["avgErr"] - lg["avgErr"])
            status = "exact" if abs(diff) < 1e-10 else "MISMATCH"
            status += " (" + how + ")"
            if not status.startswith("exact"):
                detail.append("| {} | {} | {} | {} | {:.4f} | {:.4f} | {:+.4f} | {} |".format(
                    src, rep, v, B, 100 * r["avgErr"], 100 * lg["avgErr"], diff, lg["file"]))
        counts[status] = counts.get(status, 0) + 1
    M += ["| outcome | merges |", "|---|---|"] + ["| {} | {} |".format(k, n) for k, n in sorted(counts.items())] + [""]
    M += ["Rows that differ:", "", "| source | replicate | variant | B | rebuilt % | logged % | Δ | logged file |",
          "|---|---|---|---|---|---|---|---|"] + (detail or ["| none | | | | | | | |"]) + [""]
    nolog = sorted("{}/{} {}@B{}".format(*k) for k in rebuilt if logged_match(logged, *k)[0] is None)
    M += ["Rebuilt merges with no logged row for that draw ({}): {}".format(len(nolog), ", ".join(nolog) or "none"), ""]

    # trees vs logged trees.jsonl
    M += ["## 2. Trees vs the logged trees.jsonl rows (same draw)", "",
          "| source | replicate | method | rebuilt nRF % | logged nRF % | outcome |", "|---|---|---|---|---|---|"]
    lt = {}
    for f in glob.glob(os.path.join(a.data, "logged", "*", "results_*", "trees.jsonl")):
        src = os.path.relpath(f, os.path.join(a.data, "logged")).split(os.sep)[0]
        for line in open(f):
            r = json.loads(line)
            lt[(src, r["dataset"], r["method"])] = 100 * r["RF"]
    for (src, rep, m), v in sorted(tree_rows.items()):
        logged_m = {"vote-hard-bb": "vote_hard-bb"}.get(m, m)
        x = lt.get((src, rep, logged_m))
        out = "no logged row" if x is None else ("exact" if abs(x - v) < 1e-9 else "MISMATCH")
        M.append("| {} | {} | {} | {:.2f} | {} | {} |".format(src, rep, m, v, "–" if x is None else "{:.2f}".format(x), out))
    M.append("")

    # logged summary tables recomputed from logged rows (main gcmvote draws are not banked)
    M += ["## 3. gcmvote/results/tables.md recomputed from the logged raw rows", "",
          "The REPORT's training / held-out tables use the gcmvote main branch's own MAGUS draws (only SIMMOD_R2 "
          "is from helper h5). Those draws were never banked, so they cannot be rebuilt by merging; this section "
          "recomputes the summary cells from the logged rows (data/logged/gcmvote/results/raw) to check the "
          "aggregation, and compares them cell by cell with the logged tables.md.", ""]
    main_rows = {}
    for r in logged:
        if r["impl"] == "vote" and (r["src"] == "gcmvote" or (r["src"] == "gcmvote-h5" and r["rep"] == "SIMMOD_R2"
                                                               and "vote.results" in r["file"])):
            main_rows.setdefault((r["rep"], r["variant"], r["B"]), r)
    logged_md = os.path.join(a.data, "logged", "gcmvote", "results", "tables.md")
    sections = parse_md_tables(logged_md) if os.path.exists(logged_md) else {}

    def merr(rep, v):
        r = main_rows.get((rep, v, 10))
        return None if r is None else 100 * r["avgErr"]
    nbad = ncell = 0
    for title, reps in (("Training (B = 10)", TRAIN), ("Held-out (B = 10)", HELD)):
        lt_tab = sections.get(title, {})
        for v in ["es4", "hard", "hard-bb"]:
            for lab, sel in (("proteins", lambda d: d.startswith(("BBA", "SIM"))),
                             ("DNA/RNA", lambda d: not d.startswith(("BBA", "SIM"))), ("all", lambda d: True)):
                ds = [merr(d, v) - merr(d, "magus") for d in reps if sel(d) and merr(d, v) is not None
                      and merr(d, "magus") is not None]
                mine = summary(ds, 0.05)
                theirs = lt_tab.get(("**{}**".format(lab), v))
                ncell += 1
                ok = theirs == mine
                nbad += not ok
                M.append("- {} / {} / {}: recomputed `{}`, logged `{}` {}".format(title, lab, v, mine, theirs,
                                                                                    "" if ok else "**MISMATCH**"))
    M += ["", "{} of {} summary cells differ.".format(nbad, ncell), ""]

    open(os.path.join(a.out, "TABLES.md"), "w").write("\n".join(T) + "\n")
    open(os.path.join(a.out, "MISMATCHES.md"), "w").write("\n".join(M) + "\n")
    print("wrote", os.path.join(a.out, "TABLES.md"), "and MISMATCHES.md;", counts)


def parse_md_tables(path):
    """{section title: {(row label, column name): cell}} for the '## ...' sections of a generated tables.md."""
    out, title, header = {}, None, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("## "):
            title, header = line[3:].strip(), None
        elif line.startswith("|") and title:
            cells = [c.strip() for c in line.strip("|").split("|")]
            if header is None:
                header = cells
            elif not set(line) <= set("|-"):
                for h, c in zip(header[2:], cells[2:]):
                    out.setdefault(title, {})[(cells[0], h)] = c
    return out


if __name__ == "__main__":
    main()
