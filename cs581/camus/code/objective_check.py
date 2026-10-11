"""Is the CAMUS objective aligned with accuracy? For each replicate, compare the quartet
objective (filtered gene-tree quartets displayed, t=0.5) of: the true network restricted to
its best 1-reticulation subnetwork proxy (all displayed trees of the true network), the
default CAMUS k=1 network, and the k=1 network built on the true major tree.

usage: python3 objective_check.py RESDIR COND GT OUT_TSV
"""
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import netmetric as nm  # noqa: E402
import quartets as qq  # noqa: E402

DATA = "/opt/data/camus/sim/camus-dataset"
resdir, cond, gt, out = sys.argv[1:5]
with open(out, "w") as fo:
    fo.write("rep\tr_true\tscore_true_net\tscore_default\tscore_true_major_k1\terr_default\terr_true_major\n")
    for f in sorted(glob.glob(f"{resdir}/{cond}_{gt}_*.jsonl")):
        recs = {json.loads(l)["variant"]: json.loads(l) for l in open(f)}
        if "default" not in recs or "true_major" not in recs:
            continue
        rep = recs["default"]["rep"]
        tn = open(f"{DATA}/{cond}/{rep}/true_net.nwk").readline()
        troot, _ = nm.parse_enewick(tn)
        qt = qq.QuartetTable(sorted(nm.leaves(troot)))
        C = qt.counts([l for l in open(f"{DATA}/{cond}/{rep}/{gt}.nwk") if l.strip()])
        M = qq.filter_mask(C, 0.5)
        s_true = qq.network_score(qt, C, M, nm.displayed_tree_newicks(troot))
        e = lambda r: (r["fn_k1"] + r["fp_k1"]) / 2  # noqa: E731
        fo.write(f"{rep}\t{recs['default']['r_true']}\t{s_true}\t{recs['default']['score_k1']}\t"
                 f"{recs['true_major']['score_k1']}\t{e(recs['default']):.4f}\t{e(recs['true_major']):.4f}\n")
        fo.flush()
