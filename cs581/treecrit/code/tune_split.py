# AI-assisted (Claude), exploration code for CS581 project
"""Alignment-only evaluation of splitcols.py thresholds on training draws (SIMHIGH R1-R10).

    python3 tune_split.py BASE TAU,TAU,... KEY [...]   -> data/tune_split.jsonl

Cost = Part A within-draw price model (estimated + true alignment rows, SIMHIGH):
    predicted RF change = 0.028 * dSPFN + 0.151 * dSPFP   (points; vs the unsplit BASE)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from measures import measures, W  # noqa: E402
from splitcols import split  # noqa: E402

PFN, PFP = 0.028, 0.151
base, taus, keys = sys.argv[1], [float(t) for t in sys.argv[2].split(",")], sys.argv[3:]
with open(os.path.join(HERE, "..", "data", "tune_split.jsonl"), "a") as f:
    for key in keys:
        vd = os.path.join(W, key, "vote", base)
        b = measures(key, base, os.path.join(vd, "out.fasta"))
        for tau in taus:
            out, L, Ln = split(vd, tau, os.path.join(vd, "out.split{}.fasta".format(tau)))
            m = measures(key, "{}+split{}".format(base, tau), out)
            row = {"key": key, "base": base, "tau": tau, "cols": L, "cols_split": Ln,
                   "dSPFN": 100 * (m["SPFN"] - b["SPFN"]), "dSPFP": 100 * (m["SPFP"] - b["SPFP"]),
                   "dsplit_res": 100 * (m["split_res"] - b["split_res"]),
                   "dmisplaced_res": 100 * (m["misplaced_res"] - b["misplaced_res"])}
            row["pred_dRF"] = PFN * row["dSPFN"] + PFP * row["dSPFP"]
            f.write(json.dumps(row) + "\n"); f.flush()
            print(key, tau, Ln - L, "dSPFN %+.2f dSPFP %+.2f pred dRF %+.3f" % (row["dSPFN"], row["dSPFP"], row["pred_dRF"]), flush=True)
