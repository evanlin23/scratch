"""Identifiability-aware output for a DecoDiPhy result.
usage: python decodiphy_idcheck.py LABELED_TREE.nwk all_rounds.json [round_index (default: final)]
Prints the flags (boundary / adjacent / Laplacian kernel), the per-placement abundance interval
[lo, hi] over all exact-equivalent solutions on the fitted edges, and a canonical representative."""
import sys, json
from idclass import URTree
R = URTree(open(sys.argv[1]).read().strip().split("\n")[0])
rounds = json.load(open(sys.argv[2]))
r = rounds[int(sys.argv[3])] if len(sys.argv) > 3 else [x for x in rounds if x["rounds"] == "final"][-1]
pls = [R.placement(a, x) for a, x in zip(r["anchors"], r["x"])]
fl = R.flags(pls, r["x"]); fl.pop("W")
C = R.eq_class(pls, r["p"], r["y"])
ident = not (fl["adjacent"] or fl["kernel"] or fl["boundary"]) and (C["hi"] - C["lo"]).max() < 1e-6
print("identifiable:", ident, fl)
print("anchor\tp(DecoDiPhy)\tp_lo\tp_hi\tp_canonical")
for a, p, lo, hi, pc in zip(r["anchors"], r["p"], C["lo"], C["hi"], C["canon_p"]):
    print(f"{a}\t{p:.4f}\t{lo:.4f}\t{hi:.4f}\t{pc:.4f}")
print(f"ybar: DecoDiPhy {r['y']:.4f}, canonical {C['canon_y']:.4f}")
