# AI-assisted (Claude), exploration code for CS581 project
"""Pairwise-distance and close-relative measures (added after the first Part A pass).

    python3 measures2.py [--lanes 4]  -> data/measures2.jsonl (restartable)

Taxon pairs: every cherry of the true tree, plus 3,000 random pairs (seed 1, same pairs for every alignment
of a draw). For a pair (a, b) and an alignment: compared sites = columns where both have a residue;
p = fraction of compared sites that differ.
  dist_mae_rand / dist_mae_cherry   mean |p_est - p_true| in % points (random pairs / cherries)
  dist_bias_rand / dist_bias_cherry mean (p_est - p_true)
  dist_rank_err                     1 - Spearman(p_est, p_true) over random pairs
  sites_rand                        compared sites est / true (random pairs)
  spfp_cherry / spfn_cherry         SP error restricted to residue pairs between the two taxa of a cherry
  spfp_rand / spfn_rand             the same for random pairs
"""
import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from measures import Tree, read_fasta, tree_file, W  # noqa: E402

DATA = os.path.join(HERE, "..", "data")


def rows_of(aln, taxa):
    return np.array([np.frombuffer(aln[t].upper().encode(), dtype=np.uint8) for t in taxa])


def rescol(R):
    """per row: column index of each residue."""
    g = (R == ord("-")) | (R == ord("."))
    return [np.nonzero(~g[i])[0] for i in range(R.shape[0])]


def pairs_for(key, taxa):
    tree = Tree(tree_file(os.path.join(W, key)), taxa)
    idx = {t: i for i, t in enumerate(taxa)}
    cher = []
    for v, ch in enumerate(tree.children):
        leaves = [c for c in ch if not tree.children[c]]
        if len(ch) == 2 and len(leaves) == 2:
            cher.append((idx[tree.names[leaves[0]]], idx[tree.names[leaves[1]]]))
    rng = np.random.default_rng(1)
    a = rng.integers(0, len(taxa), 3000); b = rng.integers(0, len(taxa) - 1, 3000); b = b + (b >= a)
    return cher, list(zip(a, b))


def pair_stats(R, C, pairs, Rt, Ct):
    """p-distances and residue-pair SP counts for each taxon pair."""
    gapR = (R == ord("-")) | (R == ord("."))
    gapT = (Rt == ord("-")) | (Rt == ord("."))
    out = []
    for a, b in pairs:
        both = ~gapR[a] & ~gapR[b]
        n = both.sum()
        p = (R[a][both] != R[b][both]).mean() if n else np.nan
        bt = ~gapT[a] & ~gapT[b]
        pt = (Rt[a][bt] != Rt[b][bt]).mean()
        # residue pairs (i of a, j of b) aligned in est / true
        # map est column -> residue index for a and b
        ea = np.full(R.shape[1], -1); ea[C[a]] = np.arange(len(C[a]))
        eb = np.full(R.shape[1], -1); eb[C[b]] = np.arange(len(C[b]))
        m = (ea >= 0) & (eb >= 0)
        est_pairs = set(zip(ea[m].tolist(), eb[m].tolist()))
        ta = np.full(Rt.shape[1], -1); ta[Ct[a]] = np.arange(len(Ct[a]))
        tb = np.full(Rt.shape[1], -1); tb[Ct[b]] = np.arange(len(Ct[b]))
        m = (ta >= 0) & (tb >= 0)
        true_pairs = set(zip(ta[m].tolist(), tb[m].tolist()))
        sh = len(est_pairs & true_pairs)
        out.append((p, pt, n, bt.sum(), sh, len(est_pairs), len(true_pairs)))
    return np.array(out, dtype=float)


def run(arg):
    key, method, path = arg
    try:
        true = read_fasta(os.path.join(W, key, "true.fasta"))
        taxa = sorted(true)
        est = read_fasta(path)
        Rt, R = rows_of(true, taxa), rows_of(est, taxa)
        Ct, C = rescol(Rt), rescol(R)
        if path.endswith(".masked.fasta"):
            # masked rows lose residues; residue indices no longer match -> SP parts use the unmasked file
            full = read_fasta(path.replace(".masked.fasta", ".fasta"))
            Rf = rows_of(full, taxa); Cf = rescol(Rf)
        else:
            Rf, Cf = R, C
        cher, rnd = pairs_for(key, taxa)
        r = {"key": key, "method": method}
        for lab, pairs in (("cherry", cher), ("rand", rnd)):
            s = pair_stats(R, C, pairs, Rt, Ct)
            sp = pair_stats(Rf, Cf, pairs, Rt, Ct) if Rf is not R else s
            r["dist_mae_" + lab] = 100 * float(np.nanmean(np.abs(s[:, 0] - s[:, 1])))
            r["dist_bias_" + lab] = 100 * float(np.nanmean(s[:, 0] - s[:, 1]))
            r["spfp_" + lab] = 100 * (1 - sp[:, 4].sum() / max(sp[:, 5].sum(), 1))
            r["spfn_" + lab] = 100 * (1 - sp[:, 4].sum() / sp[:, 6].sum())
            if lab == "rand":
                r["dist_rank_err"] = float(1 - stats.spearmanr(s[:, 0], s[:, 1], nan_policy="omit")[0])
                r["sites_rand"] = float(s[:, 2].sum() / s[:, 3].sum())
        return r
    except Exception as e:  # noqa: BLE001
        import traceback
        return {"key": key, "method": method, "error": repr(e), "tb": traceback.format_exc()[-1200:]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lanes", type=int, default=4)
    ap.add_argument("--paths", default=os.path.join(DATA, "aln_paths.jsonl"))
    ap.add_argument("--out", default=os.path.join(DATA, "measures2.jsonl"))
    a = ap.parse_args()
    jobs = [(j["key"], j["method"], j["path"]) for j in map(json.loads, open(a.paths))]
    keys = sorted({k for k, _, _ in jobs})
    jobs += [(k, "true", os.path.join(W, k, "true.fasta")) for k in keys]
    done = set()
    if os.path.exists(a.out):
        done = {(r["key"], r["method"]) for r in map(json.loads, open(a.out)) if "error" not in r}
    jobs = [j for j in jobs if (j[0], j[1]) not in done]
    print(len(jobs), "jobs", flush=True)
    with Pool(a.lanes) as p, open(a.out, "a") as f:
        for r in p.imap_unordered(run, jobs):
            if "error" in r:
                print("ERR", r["key"], r["method"], r["error"], r["tb"], flush=True)
                continue
            f.write(json.dumps(r) + "\n"); f.flush()


if __name__ == "__main__":
    main()
