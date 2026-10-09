"""How much can ANY merger that keeps every subset tree T_i induced improve on GTM?

For each true split B missing from the GTM tree, B is "recoverable" only if its
restriction to every S_i is trivial or a split of T_i. Splits that conflict with a
subset tree can never be recovered by a disjoint tree merger (blended or not).
Also report subset-tree error (true tree restricted to S_i vs T_i).
Usage: python3 ceiling.py OUT_TSV
"""
import sys
from datasets import REPS, inputs
from phylo import read_tree, bipartitions, fn_fp


def restrict_norm(b, mask):
    r = b & mask
    c = bin(r).count("1")
    m = bin(mask).count("1")
    if c <= 1 or c >= m - 1:
        return None
    low = mask & -mask
    return r if not (r & low) else mask ^ r


def main():
    rows = []
    for cond, reps in REPS.items():
        for rep in reps:
            for g in ("FT", "IQ"):
                x = inputs(cond, rep, g)
                T = read_tree(x["true"])
                G = read_tree(x["published_gtm"])
                names = sorted(T.label.values())
                idx = {n: i for i, n in enumerate(names)}
                bt, bg = bipartitions(T, idx), bipartitions(G, idx)
                subs = []
                for p in x["subsets"]:
                    ti = read_tree(p)
                    mask = 0
                    for v in ti.leaves():
                        mask |= 1 << idx[ti.label[v]]
                    # T_i splits as global bitsets normalized within mask
                    sp = {restrict_norm(b, mask) for b in bipartitions(ti, idx)}
                    subs.append((mask, sp))
                missing = bt - bg
                rec = 0
                for b in missing:
                    if all((r := restrict_norm(b, m)) is None or r in sp for m, sp in subs):
                        rec += 1
                # true splits incompatible with some T_i (lost by any DTM)
                lost = sum(1 for b in bt if not all((r := restrict_norm(b, m)) is None or r in sp
                                                    for m, sp in subs))
                # subset-tree error
                sub_fn = 0
                sub_tot = 0
                for p in x["subsets"]:
                    ti = read_tree(p)
                    f = fn_fp(ti, T)
                    sub_fn += f[2]
                    sub_tot += f[4]
                n = len(bt)
                rows.append((cond, rep, g, len(missing) / n, (len(missing) - rec) / n, lost / n, rec,
                             sub_fn / sub_tot))
                print(*rows[-1], sep="\t", flush=True)
    with open(sys.argv[1], "w") as f:
        f.write("condition\treplicate\tguide\tFN_GTM\tFN_floor_any_DTM_given_GTM\tFN_floor_any_DTM\t"
                "n_recoverable_missing\tsubset_tree_FN\n")
        for r in rows:
            f.write("%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%d\t%.5f\n" % r)
    import collections
    agg = collections.defaultdict(list)
    for r in rows:
        agg[(r[0], r[2])].append(r)
    print("\ncond guide FN_GTM  floor(any DTM)  recoverable/rep  subsetFN")
    for k, v in agg.items():
        print(k, "%.2f%%  %.2f%%  %.1f  %.2f%%" % (100 * sum(r[3] for r in v) / len(v),
                                                   100 * sum(r[5] for r in v) / len(v),
                                                   sum(r[6] for r in v) / len(v),
                                                   100 * sum(r[7] for r in v) / len(v)))


main()
