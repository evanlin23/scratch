"""Tables, time model, pre-registered selection and plot from results/state/*.json.

    python analyze.py > ../results/summary.md
"""
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
sys.path.insert(0, HERE)
from mf_bench import K_of, parse_setting, setting_name  # noqa: E402

TRAIN = ["BBA0101_R0", "BBA0190_R0", "SIMHIGH_R1", "1000M2_R0", "1000L1_R0", "RNASim1000_R0"]
HELD = ["BBA0067_R0", "SIMHIGH_R2", "1000M3_R0", "RNASim1000_R1"]
SETTINGS = [setting_name(n, s, p) for s in (200, 100) for n in (10, 6, 4, 3) for p in (0, 1)]
TIE = 0.05


def load():
    return {os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob(os.path.join(RES, "state", "*.json"))}


def pct(x, tot):
    return "{:.0f} ({:.0f}%)".format(x, 100 * x / tot) if tot else "-"


def profile_table(S):
    print("## Q1. Where MAGUS spends its time (4 cores, pure MAGUS, paper flags)\n")
    print("Seconds (share of end-to-end wall). Pool = the shared 4-slot MAFFT task pool (25 subset alignments + "
          "10 backbones), split by measured MAFFT CPU share. Graph build / MCL / trace+write = the time after the "
          "last MAFFT task finished.\n")
    print("| dataset | type | wall s | guide tree + decomposition | subset alignments | backbones | graph build | MCL | trace + write | backbone share of MAFFT CPU |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for name in TRAIN + HELD:
        st = S.get(name, {})
        if "default" not in st:
            continue
        p = st["default"]["prof"]
        tot = p["total"]
        cs, cb = p["subset_cpu"], p["backbone_cpu"]
        pool = p["pool_end"] - p["decomposition"]
        sub, bb = pool * cs / (cs + cb), pool * cb / (cs + cb)
        graph = p["graph_exposed"] + p["graph_write"]
        rest = tot - p["pool_end"] - graph - p["mcl"]
        print("| {} | {} | {:.0f} | {} | {} | {} | {} | {} | {} | {:.0f}% |".format(
            name + (" (held out)" if name in HELD else ""), st["datatype"], tot, pct(p["decomposition"], tot),
            pct(sub, tot), pct(bb, tot), pct(graph, tot), pct(p["mcl"], tot), pct(rest, tot), 100 * cb / (cs + cb)))
    print()
    print("### Graph build: MAGUS's Python builder vs vectorized (merge-only, same 25 subsets and 10 backbones)\n")
    print("| dataset | graph entries | identical graph | identical clusters / output | Python build s | vectorized build s | MCL 1 thread s | MCL 4 threads s | identical with MCL -te 4 | merge wall: pure MAGUS -> both exact speed-ups |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for name in TRAIN + HELD:
        g = S.get(name, {}).get("graph")
        if not g:
            continue
        py, fg, fm = g["python"], g["fast"], g["fast_mcl4"]
        print("| {} | {:,} | {} | {} / {} | {:.1f} | {:.1f} | {:.1f} | {:.1f} | {} | {:.0f} -> {:.0f} s |".format(
            name, g["graph_entries"], g["graph_identical"], g["clusters_identical_fast"], g["output_identical_fast"],
            py["prof"]["graph_buildMatrix"] + py["prof"]["graph_write"], fg["prof"]["graph_buildMatrix"] + fg["prof"]["graph_write"],
            fg["prof"]["mcl"], fm["prof"]["mcl"], g["output_identical_mcl4"], py["merge_wall"], fm["merge_wall"]))
    print()


def predicted(st, setting):
    """Predicted end-to-end wall of a setting (seconds) from the default run's profile."""
    p = st["default"]["prof"]
    n, s, pr = parse_setting(setting)
    cs, cb_each = p["subset_cpu"], p["bb_cpu_each"]
    rate = (cs + sum(cb_each)) / (p["pool_end"] - p["decomposition"])  # CPU-s per wall-s in the pool
    if s == 200:
        cb = sum(cb_each[:n])
    else:
        cb = sum(st["bb100"]["cpu_each"][:n])
    merge = st["sweep"][setting]["merge_wall"] if setting in st["sweep"] else st["graph"]["fast_mcl4"]["merge_wall"]
    return p["decomposition"] + (cs + cb) / rate + merge


def predicted_default(st):
    p = st["default"]["prof"]
    return p["pool_end"] + st["graph"]["python"]["merge_wall"]


def err(x):
    return x["err"]


def sweep_tables(S):
    print("## Q2. Fewer / smaller backbones, with and without support pruning (merge-only, paired)\n")
    print("Δ = setting SP error − default merge (10 x 200, no pruning) on the same subsets, in points; "
          "time ratio = predicted end-to-end wall / measured default MAGUS wall. Pruning keeps edges with "
          "support >= K = ceil(0.4 N).\n")
    names = [n for n in TRAIN if "sweep" in S.get(n, {})]
    print("| setting | K | " + " | ".join(names) + " | mean Δ | W/T/L | max Δ | mean ΔSPFN | mean ΔSPFP | time ratio (geo-mean) |")
    print("|---|---|" + "---|" * len(names) + "---|---|---|---|---|---|")
    rows = {}
    for setting in SETTINGS:
        n, s, p = parse_setting(setting)
        ds, dfn, dfp, tr = [], [], [], []
        for name in names:
            st = S[name]
            base = st["graph"]["fast_mcl4"]
            x = st["sweep"].get(setting) if setting != "n10s200p0" else base
            if x is None:
                ds.append(None)
                continue
            ds.append(err(x) - err(base))
            dfn.append(100 * (x["SPFN"] - base["SPFN"]))
            dfp.append(100 * (x["SPFP"] - base["SPFP"]))
            tr.append(predicted(st, setting) / st["default"]["wall"])
        v = [d for d in ds if d is not None]
        if not v:
            continue
        w = sum(d < -TIE for d in v)
        l = sum(d > TIE for d in v)
        geo = math.exp(sum(map(math.log, tr)) / len(tr))
        rows[setting] = {"d": ds, "max": max(v), "mean": sum(v) / len(v), "time": geo, "complete": len(v) == len(names)}
        print("| {} | {} | {} | {:+.2f} | {}/{}/{} | {:+.2f} | {:+.2f} | {:+.2f} | {:.2f} ({:.2f}x) |".format(
            setting, K_of(n) if p else "-", " | ".join("-" if d is None else "{:+.2f}".format(d) for d in ds),
            sum(v) / len(v), w, len(v) - w - l, l, max(v), sum(dfn) / len(dfn), sum(dfp) / len(dfp), geo, 1 / geo))
    print()
    return rows, names


def main():
    S = load()
    profile_table(S)
    if not any("sweep" in S.get(n, {}) for n in TRAIN):
        return
    rows, names = sweep_tables(S)
    print("### Time model check (default): predicted vs measured wall\n")
    print("| dataset | measured s | predicted s |\n|---|---|---|")
    for name in names:
        print("| {} | {:.0f} | {:.0f} |".format(name, S[name]["default"]["wall"], predicted_default(S[name])))
    print()
    complete = len(names) == len(TRAIN)
    elig = [s for s, r in rows.items() if r["complete"] and r["max"] <= 0.25]
    print("## Q3. Pre-registered choice{}\n".format("" if complete else " (PROVISIONAL: training incomplete)"))
    print("Eligible (Δ <= +0.25 on every training dataset): " + ", ".join(
        "{} ({:.2f}x)".format(s, 1 / rows[s]["time"]) for s in sorted(elig, key=lambda s: rows[s]["time"])))
    if elig:
        best = min(elig, key=lambda s: (round(rows[s]["time"] / 0.02), -parse_setting(s)[0]))
        print("\n**Chosen: {}** (predicted {:.2f}x faster).".format(best, 1 / rows[best]["time"]))
    sec = [s for s, r in rows.items() if r["complete"] and r["max"] <= 0.5]
    if sec:
        b2 = min(sec, key=lambda s: rows[s]["time"])
        print("\nSecondary (post hoc): fastest within +0.5 everywhere: {} ({:.2f}x).".format(b2, 1 / rows[b2]["time"]))
    print()
    json.dump(rows, open(os.path.join(RES, "sweep_rows.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
