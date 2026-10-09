"""Adaptive soft constraints: relax only the subset alignments the evidence disputes.

    python -m gcmx.adaptive_split SUBALIGNMENT_DIR EVIDENCE_DIR OUT_DIR
        [--rule tiers|tiers-hard|proportional] [--report REPORT.json]

Disagreement of subset alignment S with the evidence (no reference needed):
over the evidence alignments (e.g. HMM-extended backbones, which contain
every sequence), the fraction of letter pairs of S that the evidence puts in
one column but S puts in different columns. On 1000M2 R0 this predicts the
true SP error of S with Spearman rho = 0.68.

Rules (q = quantiles of disagreement over the subsets of this replicate):
  tiers        d < q40: 1 group; q40-q80: 3 groups; > q80: 6 groups
  tiers-hard   d < q60: 1 group; q60-q90: 3 groups; > q90: 1 group per ~4 sequences
  proportional groups = clip(round(3 * d / median(d)), 1, 8)
Groups are similarity groups (gcmx.split, average linkage on p-distance).
"""

import argparse
import collections
import json
import os
import random

import numpy as np

from . import fasta
from .compact import letter_columns
from .split import split_groups


def disagreement(sub, evidence_files):
    names = list(sub)
    subcol = {n: np.array([c for c, ch in enumerate(sub[n]) if ch not in "-."]) for n in names}
    split_pairs = ev_pairs = 0
    for path in evidence_files:
        ev = letter_columns(path, names)
        both, per_ev = collections.Counter(), collections.Counter()
        for n in names:
            if n not in ev:
                continue
            e, s = ev[n], subcol[n]
            m = e >= 0
            for sc, ec in zip(s[m].tolist(), e[m].tolist()):
                both[sc, ec] += 1
                per_ev[ec] += 1
        same_ev = sum(v * (v - 1) / 2 for v in per_ev.values())
        same_both = sum(v * (v - 1) / 2 for v in both.values())
        split_pairs += same_ev - same_both
        ev_pairs += same_ev
    return split_pairs / ev_pairs if ev_pairs else 0.0


def groups_for(d, ds, n, rule):
    q = lambda p: float(np.quantile(ds, p))
    if rule == "tiers":
        return 1 if d < q(0.4) else 3 if d < q(0.8) else 6
    if rule == "tiers-hard":
        return 1 if d < q(0.6) else 3 if d < q(0.9) else max(1, n // 4)
    return int(np.clip(round(3 * d / float(np.median(ds))), 1, 8))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("subalignments")
    parser.add_argument("evidence")
    parser.add_argument("outdir")
    parser.add_argument("--rule", default="tiers", choices=("tiers", "tiers-hard", "proportional"))
    parser.add_argument("--report", default=None)
    args = parser.parse_args()

    evidence_files = [os.path.join(args.evidence, f) for f in sorted(os.listdir(args.evidence))]
    files = sorted(os.listdir(args.subalignments))
    subs = {f: fasta.read(os.path.join(args.subalignments, f)) for f in files}
    ds = {f: disagreement(subs[f], evidence_files) for f in files}
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(0)
    plan = {}
    for f in files:
        m = groups_for(ds[f], list(ds.values()), len(subs[f]), args.rule)
        plan[f] = {"disagreement": round(ds[f], 4), "groups": m}
        for g, taxa in enumerate(split_groups(subs[f], m, "linkage", rng)):
            fasta.write(fasta.restrict(subs[f], taxa), os.path.join(args.outdir, "{}_g{}.txt".format(f[:-4], g)))
    if args.report:
        json.dump(plan, open(args.report, "w"), indent=1)
    print("groups per subset:", sorted(p["groups"] for p in plan.values()), "->", len(os.listdir(args.outdir)), "groups")


if __name__ == "__main__":
    main()
