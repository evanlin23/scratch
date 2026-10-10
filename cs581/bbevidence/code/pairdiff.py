"""Pair-level comparison of two backbone tools on the SAME sequence sets.

    python3 pairdiff.py REP_DIR TOOL_A TOOL_B   -> one JSON row appended to REP_DIR/pairdiff.jsonl

For each of the 10 backbones, the cross-subset residue pairs aligned by A and by B are split into
shared (A and B), A-only and B-only, each with its precision against the reference, plus where the
A-only false pairs sit (gappy reference columns, low-identity regions). Also a reference-free
statistic usable as a predictor: agreement = |shared| / |A|.
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bbe  # noqa: E402


def pairs_of(R, aln, gid):
    """Cross-subset aligned residue pairs as sorted int64 keys gid_x * G + gid_y (gid_x < gid_y),
    and whether each is a true pair."""
    rc = bbe.residue_columns(aln)
    cols, ids, refs, subs = [], [], [], []
    for t in aln:
        cols.append(rc[t]); ids.append(gid[t]); refs.append(R.refcol[t])
        subs.append(np.full(len(rc[t]), R.subset_of[t], dtype=np.int32))
    cols, ids, refs, subs = map(np.concatenate, (cols, ids, refs, subs))
    order = np.argsort(cols, kind="stable")
    cols, ids, refs, subs = cols[order], ids[order], refs[order], subs[order]
    bounds = np.flatnonzero(np.diff(cols)) + 1
    G = int(max(v.max() for v in gid.values())) + 1
    ks, ts, gap = [], [], []
    for lo, hi in zip(np.r_[0, bounds], np.r_[bounds, len(cols)]):
        if hi - lo < 2:
            continue
        i, j = np.triu_indices(hi - lo, 1)
        i, j = i + lo, j + lo
        keep = subs[i] != subs[j]
        i, j = i[keep], j[keep]
        a, b = np.minimum(ids[i], ids[j]), np.maximum(ids[i], ids[j])
        ks.append(a.astype(np.int64) * G + b)
        ts.append(refs[i] == refs[j])
        gap.append(np.maximum(R.ref_gapfrac[refs[i]], R.ref_gapfrac[refs[j]]))
    k, t, g = np.concatenate(ks), np.concatenate(ts), np.concatenate(gap)
    o = np.argsort(k)
    return k[o], t[o], g[o]


def main(rep, ta, tb):
    R = bbe.Rep(rep)
    gid, off = {}, 0
    for t, rc in R.refcol.items():
        gid[t] = np.arange(off, off + len(rc), dtype=np.int64)
        off += len(rc)
    A, _, _ = bbe.parse_variant(rep, ta)
    B, _, _ = bbe.parse_variant(rep, tb)
    acc = {k: 0 for k in ("shared", "shared_tp", "a_only", "a_only_tp", "b_only", "b_only_tp",
                          "a_only_fp_gappy", "b_only_fp_gappy", "a", "b", "a_only_gappy", "b_only_gappy")}
    for (_, a), (_, b) in zip(A, B):
        assert set(a) == set(b)
        ka, tA, gA = pairs_of(R, a, gid)
        kb, tB, gB = pairs_of(R, b, gid)
        inb = np.isin(ka, kb, assume_unique=True)
        ina = np.isin(kb, ka, assume_unique=True)
        acc["a"] += len(ka); acc["b"] += len(kb)
        acc["shared"] += int(inb.sum()); acc["shared_tp"] += int(tA[inb].sum())
        acc["a_only"] += int((~inb).sum()); acc["a_only_tp"] += int(tA[~inb].sum())
        acc["b_only"] += int((~ina).sum()); acc["b_only_tp"] += int(tB[~ina].sum())
        acc["a_only_fp_gappy"] += int(((~inb) & ~tA & (gA > 0.5)).sum())
        acc["b_only_fp_gappy"] += int(((~ina) & ~tB & (gB > 0.5)).sum())
        acc["a_only_gappy"] += int(((~inb) & (gA > 0.5)).sum())
        acc["b_only_gappy"] += int(((~ina) & (gB > 0.5)).sum())
    row = {"rep": os.path.basename(rep.rstrip("/")), "A": ta, "B": tb,
           "agreement": round(acc["shared"] / acc["a"], 4),
           "prec_shared": round(acc["shared_tp"] / max(acc["shared"], 1), 4),
           "prec_a_only": round(acc["a_only_tp"] / max(acc["a_only"], 1), 4),
           "prec_b_only": round(acc["b_only_tp"] / max(acc["b_only"], 1), 4),
           "frac_a_only": round(acc["a_only"] / acc["a"], 4),
           "frac_b_only": round(acc["b_only"] / acc["b"], 4),
           "a_only_fp_share_gappy": round(acc["a_only_fp_gappy"] / max(acc["a_only"] - acc["a_only_tp"], 1), 4),
           "b_only_fp_share_gappy": round(acc["b_only_fp_gappy"] / max(acc["b_only"] - acc["b_only_tp"], 1), 4),
           "a_only_share_gappy": round(acc["a_only_gappy"] / max(acc["a_only"], 1), 4),
           "b_only_share_gappy": round(acc["b_only_gappy"] / max(acc["b_only"], 1), 4),
           "n_a": acc["a"], "n_b": acc["b"]}
    with open(os.path.join(rep, "pairdiff.jsonl"), "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row))


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), sys.argv[2], sys.argv[3])
