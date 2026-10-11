"""Exact-equivalence test of forest_fast vs distphylo's own functions on small simulated data
(n = 24-30), over many (m, M, tau) points including small M (ball smaller than the component, so
the Extender is exercised). Compares the full split sets per component and the validity flag.

usage: python validate_small.py DISTPHYLO_DIR OUT_TSV
"""
import sys, os, csv, tempfile, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from validate_forest import distphylo_forest, nontrivial_from_newick  # noqa: E402 (adds DP to path)
import forest_fast as F  # noqa: E402
import pipeline as P  # noqa: E402
from common import read_fasta, jc_distance  # noqa: E402


def main():
    out = sys.argv[2]
    rows = []
    w = tempfile.mkdtemp()
    for seed, (n, k, reg) in enumerate([(24, 300, "U:0.05:0.1"), (30, 200, "U:0.005:0.05"),
                                         (30, 300, "U:0.1:0.3"), (26, 500, "UH:1.0:0.7")]):
        tf, aln = P.simulate(n, k, reg, "JC", 100 + seed, w)
        names, X = read_fasta(aln)
        D = jc_distance(X)
        fin = D[np.isfinite(D)]
        b = P.mst_bottleneck(D)
        for mf, Mf, tq in itertools.product((1.0001, 0.7, 0.45), (2.05, 3.0, 6.0), (0.05, 0.2, 0.5)):
            m = b * mf
            tau = float(np.quantile(fin[fin > 0], 0.01)) * tq
            if not (m > 3 * tau):
                continue
            M = max(2 * m + 3 * tau + 1e-3, m * Mf / 2.0 * 1.0) if Mf < 3 else m * Mf
            if not (M > 2 * m + 3 * tau):
                M = 2 * m + 3 * tau + 1e-3
            dp = distphylo_forest(D, m, M, tau)
            ff = F.run_forest(D, m, M, [tau])[tau]
            if dp is None:
                dp_splits = None
            else:
                dp_splits = sorted(sorted(map(int, s)) for s in
                                   [nontrivial_from_newick(nw, c) if nw else set() for c, nw in dp])
            ff_splits = sorted(sorted(map(int, s)) for s in ff["splits"]) if ff["valid"] else None
            ball_full = bool(np.all(D[np.isfinite(D)] < M))
            rows.append(dict(seed=seed, n=n, k=k, regime=reg, m=round(m, 4), M=round(M, 4), tau=round(tau, 5),
                             ncomp=len(ff["comps"]), ball_covers_all=ball_full,
                             distphylo_valid=dp is not None, fast_valid=ff["valid"],
                             distphylo_nsplits=None if dp is None else sum(len(s) for s in dp_splits),
                             fast_nsplits=sum(len(s) for s in ff["splits"]),
                             identical=(dp_splits == ff_splits)))
            print(rows[-1], flush=True)
    with open(out, "w") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)
    print("identical:", sum(r["identical"] for r in rows), "/", len(rows))


if __name__ == "__main__":
    main()
