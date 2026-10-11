"""FN / FP / RF rates of estimated trees vs a reference tree (unrooted bipartitions).

    python treeerr.py TRUE_TREE EST_TREE [EST_TREE ...]

Pure Python (bitmask bipartitions), restricted to the common leaf set. FN rate =
missing true internal bipartitions / true internal bipartitions; FP rate likewise
over the estimated tree. Matches cs581/ml/code/treeerr.py (dendropy) on binary trees.
"""
import re
import sys

TOKEN = re.compile(r"\(|\)|,|;|:[^,();]*|\[[^\]]*\]|[^(),;:\[\]]+")


def leaf_sets(newick):
    """Return (labels, list of leaf-label sets below each internal node)."""
    stack, clusters, labels, last = [[]], [], [], None
    for tok in TOKEN.findall(newick):
        tok = tok.strip()
        if not tok or tok[0] in ":[":
            continue
        if tok == "(":
            stack.append([])
        elif tok == ")":
            kids = stack.pop()
            clusters.append(kids)
            stack[-1].extend(kids)
            last = ")"
        elif tok in ",;":
            last = tok
        else:
            if last == ")":  # internal node label / support value
                last = None
                continue
            labels.append(tok.strip("'\""))
            stack[-1].append(tok.strip("'\""))
            last = None
    return labels, clusters


def bipartitions(newick, keep):
    labels, clusters = leaf_sets(newick)
    order = {x: i for i, x in enumerate(sorted(keep))}
    full = (1 << len(order)) - 1
    out = set()
    for c in clusters:
        m = 0
        for x in c:
            if x in order:
                m |= 1 << order[x]
        if m & 1:
            m = full ^ m  # canonical side: the one without leaf 0
        k = bin(m).count("1")
        if 2 <= k <= len(order) - 2:
            out.add(m)
    return out


def error(true_path, est_path):
    t, e = open(true_path).read(), open(est_path).read()
    common = set(leaf_sets(t)[0]) & set(leaf_sets(e)[0])
    bt, be = bipartitions(t, common), bipartitions(e, common)
    n = len(common)
    fn, fp = len(bt - be), len(be - bt)
    return {"n": n, "fn": fn, "fp": fp, "true_int": len(bt), "est_int": len(be),
            "fn_rate": fn / len(bt) if bt else 0.0, "fp_rate": fp / len(be) if be else 0.0,
            "rf_rate": (fn + fp) / (2.0 * (n - 3))}


if __name__ == "__main__":
    for p in sys.argv[2:]:
        r = error(sys.argv[1], p)
        print("%s\tFN=%.4f\tFP=%.4f\tRF=%.4f\t(n=%d, true int=%d, est int=%d)" %
              (p, r["fn_rate"], r["fp_rate"], r["rf_rate"], r["n"], r["true_int"], r["est_int"]))
