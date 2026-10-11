# AI-assisted (Claude), exploration code for CS581 project
"""Which alignment error measure predicts FastTree error? (Part A)

    python3 analyze.py [--set simhigh|protein|all]  -> results/partA.md, results/partA.json

Units: estimated alignments (key, method) with a FastTree tree; `true` rows give each draw's floor and
oracle `split_*` rows are excluded from the main fits. RF = normalized RF to the true tree, in % points.
  within   fixed effects per draw: RF and measure demeaned within draw -> Pearson r, Spearman rho, R^2
  delta    paired vs MAGUS: dRF vs dMeasure over (draw, variant) -> Spearman rho, R^2, LODO-CV R^2
  absolute across all alignments: Spearman rho(RF, x) and rho(RF - RF_true, x)
  noise    RF differences between near-identical alignments (vote_hard vs its masked copy; pairs of
           variants within 0.25 SP points) -> the share of dRF variance no alignment measure can explain
CIs: cluster bootstrap over draws (2,000 resamples).
"""
import argparse
import json
import os

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
RES = os.path.join(HERE, "..", "results")

MEAS = [  # name, label, sign (+1: higher = more error)
    ("avgErr", "SP error (SPFN+SPFP)/2", 1), ("SPFN", "SPFN", 1), ("SPFP", "SPFP", 1),
    ("TCerr", "1 - TC (column recall error)", 1), ("colprec_err", "1 - column precision", 1),
    ("oversplit_cols", "over-split true columns", 1), ("overmerge_cols", "over-merged est. columns", 1),
    ("split_res", "split residues", 1), ("misplaced_res", "misplaced residues", 1), ("res_err", "residue error", 1),
    ("SPFN_pi", "SPFN, PI columns", 1), ("SPFP_pi", "SPFP, PI est. columns", 1),
    ("split_res_pi", "split residues, PI", 1), ("misplaced_res_pi", "misplaced residues, PI", 1),
    ("oversplit_cols_pi", "over-split PI columns", 1), ("overmerge_cols_pi", "over-merged PI est. columns", 1),
    ("pi_cols_ratio", "PI columns est/true", 1), ("len_ratio", "length est/true", 1),
    ("gap_open_ratio", "gap openings est/true", 1), ("bnd_err", "error at indel boundaries", 1),
    ("int_err", "error away from indels", 1), ("tax_err_p90", "per-taxon error, 90th pct", 1),
    ("fitch_excess", "Fitch excess on true tree", 1), ("fitch_excess_pi", "Fitch excess, PI cols", 1),
]
MAIN_VARIANTS = None  # all estimated methods


def load(which):
    trees = [json.loads(l) for l in open(os.path.join(DATA, "trees.jsonl"))]
    meas = {}
    for l in open(os.path.join(DATA, "measures.jsonl")):
        r = json.loads(l)
        meas[(r["key"], r["method"])] = r
    rftrue = {r["key"]: 100 * r["RF"] for r in trees if r["method"] == "true"}
    rows = []
    for t in trees:
        m = meas.get((t["key"], t["method"]))
        if m is None or t["method"] == "true":
            continue
        k = t["key"]
        if which == "simhigh" and not k.startswith("SIMHIGH"):
            continue
        if which == "protein" and not k.startswith("SIM"):
            continue
        r = dict(m)
        r["RF"] = 100 * t["RF"]
        r["RFtrue"] = rftrue.get(k, np.nan)
        r["TCerr"] = 1 - m["TC"]
        r["colprec_err"] = 1 - m["col_prec"]
        for x, _, _ in MEAS:
            if x in r and r[x] is not None and x not in ("len_ratio", "gap_open_ratio", "pi_cols_ratio"):
                r[x] = 100 * r[x]
        rows.append(r)
    return rows


def demean(rows, field, keys):
    v = np.array([r[field] for r in rows], dtype=float)
    out = v.copy()
    for k in set(keys):
        m = keys == k
        out[m] = v[m] - v[m].mean()
    return out


def within(rows, x, boot=2000, rng=None):
    est = [r for r in rows if not r["method"].startswith("split")]
    keys = np.array([r["key"] for r in est])
    y = demean(est, "RF", keys)
    xv = demean(est, x, keys)
    r = np.corrcoef(xv, y)[0, 1]
    rho = stats.spearmanr(xv, y)[0]
    ukeys = np.unique(keys)
    bs = []
    idx = {k: np.nonzero(keys == k)[0] for k in ukeys}
    for _ in range(boot):
        s = rng.choice(ukeys, len(ukeys))
        ii = np.concatenate([idx[k] for k in s])
        if np.std(xv[ii]) == 0:
            continue
        bs.append(np.corrcoef(xv[ii], y[ii])[0, 1])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    slope = np.polyfit(xv, y, 1)[0]
    return dict(n=len(est), draws=len(ukeys), r=r, r_lo=lo, r_hi=hi, R2=r * r, rho=rho, slope=slope)


def paired(rows, x, base="magus"):
    by = {}
    for r in rows:
        by.setdefault(r["key"], {})[r["method"]] = r
    d = []
    for k, ms in by.items():
        if base not in ms:
            continue
        for m, r in ms.items():
            if m == base or m.startswith("split"):
                continue
            d.append((k, m, r["RF"] - ms[base]["RF"], r[x] - ms[base][x]))
    return d


def lodo_r2(d):
    """leave-one-draw-out CV R^2 of dRF ~ a + b dx (vs predicting the training mean)."""
    keys = np.array([t[0] for t in d]); y = np.array([t[2] for t in d]); x = np.array([t[3] for t in d])
    pred, base = np.zeros_like(y), np.zeros_like(y)
    for k in np.unique(keys):
        tr, te = keys != k, keys == k
        b, a = np.polyfit(x[tr], y[tr], 1)
        pred[te] = a + b * x[te]
        base[te] = y[tr].mean()
    return 1 - ((y - pred) ** 2).sum() / ((y - base) ** 2).sum()


def lodo_r2_within(rows, xs):
    """LODO CV R^2 for within-draw demeaned RF ~ demeaned predictors (multiple regression, no intercept)."""
    est = [r for r in rows if not r["method"].startswith("split")]
    keys = np.array([r["key"] for r in est])
    y = demean(est, "RF", keys)
    X = np.column_stack([demean(est, x, keys) for x in xs])
    pred = np.zeros_like(y)
    for k in np.unique(keys):
        tr, te = keys != k, keys == k
        beta, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None)
        pred[te] = X[te] @ beta
    return 1 - ((y - pred) ** 2).sum() / (y ** 2).sum()


def noise(rows):
    by = {}
    for r in rows:
        by.setdefault(r["key"], {})[r["method"]] = r
    mask = [ms["vote_hard_mask"]["RF"] - ms["vote_hard"]["RF"] for ms in by.values()
            if "vote_hard_mask" in ms and "vote_hard" in ms]
    near = []
    for k, ms in by.items():
        ml = [m for m in ms if not m.startswith("split") and m != "vote_hard_mask"]
        for i in range(len(ml)):
            for j in range(i + 1, len(ml)):
                a, b = ms[ml[i]], ms[ml[j]]
                if abs(a["avgErr"] - b["avgErr"]) < 0.25 and abs(a["SPFN"] - b["SPFN"]) < 0.5:
                    near.append(a["RF"] - b["RF"])
    return np.array(mask), np.array(near)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="simhigh")
    a = ap.parse_args()
    rng = np.random.default_rng(7)
    rows = load(a.set)
    out = {"set": a.set, "n_rows": len(rows)}
    lines = ["## Part A tables ({}: {} estimated alignments with trees, {} draws)".format(
        a.set, len([r for r in rows if not r["method"].startswith("split")]), len({r["key"] for r in rows})), ""]
    # within-draw
    W = []
    for x, lab, _ in MEAS:
        if any(r.get(x) is None for r in rows):
            continue
        w = within(rows, x, rng=rng)
        d = paired(rows, x)
        dy = np.array([t[2] for t in d]); dx = np.array([t[3] for t in d])
        w["d_n"] = len(d)
        w["d_rho"], w["d_p"] = stats.spearmanr(dx, dy)
        w["d_r"] = np.corrcoef(dx, dy)[0, 1]
        w["d_R2_cv"] = lodo_r2(d)
        w["d_slope"] = np.polyfit(dx, dy, 1)[0]
        est = [r for r in rows if not r["method"].startswith("split")]
        w["abs_rho"] = stats.spearmanr([r[x] for r in est], [r["RF"] for r in est])[0]
        w["room_rho"] = stats.spearmanr([r[x] for r in est], [r["RF"] - r["RFtrue"] for r in est])[0]
        w["x"], w["label"] = x, lab
        W.append(w)
    W.sort(key=lambda w: -abs(w["r"]))
    lines += ["### Within-draw (fixed effects) and paired-vs-MAGUS association with FastTree nRF", "",
              "| measure | within r [95% CI] | within R² | within ρ | slope (RF pts per pt) | paired ΔRF~Δx ρ (p) | paired r | paired LODO-CV R² | across-draw ρ(RF) | ρ(RF − RF_true) |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for w in W:
        lines.append("| {label} | {r:+.2f} [{r_lo:+.2f}, {r_hi:+.2f}] | {R2:.2f} | {rho:+.2f} | {slope:+.3f} | {d_rho:+.2f} ({d_p:.1g}) | {d_r:+.2f} | {d_R2_cv:+.2f} | {abs_rho:+.2f} | {room_rho:+.2f} |".format(**w))
    out["within"] = W
    # two-predictor models (within, LODO CV)
    lines += ["", "### Within-draw LODO cross-validated R² of small models", "", "| predictors | CV R² |", "|---|---|"]
    combos = [["avgErr"], ["SPFN"], ["SPFP"], ["SPFN", "SPFP"], ["split_res", "misplaced_res"],
              ["SPFN_pi", "SPFP_pi"], ["oversplit_cols", "overmerge_cols"], ["fitch_excess"],
              ["fitch_excess_pi"], ["SPFN", "SPFP", "fitch_excess"], ["misplaced_res_pi"], ["SPFP_pi"],
              ["tax_err_p90"], ["TCerr"], ["SPFN", "SPFP", "TCerr", "len_ratio"]]
    out["cv"] = {}
    for c in combos:
        v = lodo_r2_within(rows, c)
        out["cv"]["+".join(c)] = v
        lines.append("| {} | {:+.2f} |".format(" + ".join(c), v))
    # noise
    mk, near = noise(rows)
    est = [r for r in rows if not r["method"].startswith("split")]
    keys = np.array([r["key"] for r in est])
    ywith = demean(est, "RF", keys)
    d = paired(rows, "avgErr")
    dy = np.array([t[2] for t in d])
    out["noise"] = dict(mask_n=len(mk), mask_sd=float(mk.std(ddof=1)) if len(mk) > 1 else None,
                        mask_mean_abs=float(np.abs(mk).mean()) if len(mk) else None,
                        near_n=len(near), near_sd=float(near.std(ddof=1)) if len(near) > 1 else None,
                        within_sd=float(ywith.std(ddof=1)), d_sd=float(dy.std(ddof=1)))
    nz = out["noise"]
    lines += ["", "### Noise floor", "",
              "- vote_hard vs its masked copy (3-8 of ~8,000 columns removed): n = {mask_n}, SD of ΔRF = {mask_sd:.2f}, mean |ΔRF| = {mask_mean_abs:.2f}".format(**nz),
              "- pairs of alignments of the same draw within 0.25 SP points (and 0.5 SPFN): n = {near_n}, SD of ΔRF = {near_sd:.2f}".format(**nz) if nz["near_sd"] else "- near pairs: too few",
              "- SD of within-draw RF deviations: {within_sd:.2f}; SD of ΔRF vs MAGUS: {d_sd:.2f}".format(**nz)]
    if nz["mask_sd"]:
        # pairwise noise SD s means each alignment carries noise var s^2/2
        ceil = 1 - (nz["mask_sd"] ** 2 / 2) / nz["within_sd"] ** 2
        out["noise"]["within_ceiling_R2"] = ceil
        lines.append("- implied ceiling on within-draw R² (if every alignment's RF carries independent noise of variance SD_mask²/2): {:.2f}".format(ceil))
    # per-variant means
    lines += ["", "### Per-variant means, paired vs MAGUS (Δ in points)", "",
              "| variant | n | ΔRF | ΔSP err | ΔSPFN | ΔSPFP | Δsplit res | Δmisplaced res | Δmisplaced PI | Δover-merged PI cols | ΔFitch excess (%) |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    by = {}
    for r in rows:
        by.setdefault(r["key"], {})[r["method"]] = r
    vs = sorted({r["method"] for r in rows} - {"magus"})
    out["variants"] = {}
    for v in vs:
        pr = [(ms[v], ms["magus"]) for ms in by.values() if v in ms and "magus" in ms]
        if not pr:
            continue
        f = lambda x: np.mean([a[x] - b[x] for a, b in pr])  # noqa: E731
        vals = dict(n=len(pr), dRF=f("RF"), dSP=f("avgErr"), dFN=f("SPFN"), dFP=f("SPFP"), dsplit=f("split_res"),
                    dmis=f("misplaced_res"), dmispi=f("misplaced_res_pi"), domp=f("overmerge_cols_pi"),
                    dfit=f("fitch_excess"))
        out["variants"][v] = vals
        lines.append("| {} | {n} | {dRF:+.2f} | {dSP:+.2f} | {dFN:+.2f} | {dFP:+.2f} | {dsplit:+.2f} | {dmis:+.2f} | {dmispi:+.2f} | {domp:+.2f} | {dfit:+.3f} |".format(v, **vals))
    os.makedirs(RES, exist_ok=True)
    open(os.path.join(RES, "partA_{}.md".format(a.set)), "w").write("\n".join(lines) + "\n")
    json.dump(out, open(os.path.join(RES, "partA_{}.json".format(a.set)), "w"), indent=1, default=float)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
