"""Rescore published species trees (FastMulRFS data, doi:10.13012/B2IDB-5721322_V1)
against the true species trees and compare with the authors' CSV.
Usage: python rescore_published.py /opt/data/gdl/fmrfs out.csv"""
import csv
import os
import sys

from phylo import read_trees, rf_error

root, out = sys.argv[1], sys.argv[2]
pub = {}
with open(os.path.join(root, "csvs/data-error-and-timings-ntaxa-100.csv")) as f:
    for r in csv.DictReader(f):
        pub[(r["PSIZ"], float(r["DLRT"]), int(r["REPL"]), r["SQLN"], r["NGEN"], r["MTHD"])] = float(r["SERF"]) if r["SERF"] != "NA" else None

files = {"astral": "astral-raxml-sqln-{s}-ngen-{n}.tree",
         "stag": "stag-raxml-sqln-{s}-ngen-{n}.tree",
         "mulrf": "mulrf-raxml-sqln-{s}-ngen-{n}.tree",
         "fastmulrfs-single": "fastmulrfs-raxml-sqln-{s}-ngen-{n}.tree.single"}
rows, diffs = [], []
for psiz in ["10000000", "50000000"]:
    for dl in ["0.0000000001", "0.0000000002", "0.0000000005"]:
        for rep in range(1, 11):
            d = os.path.join(root, "ntaxa-100.dlrate-%s.psize-%s" % (dl, psiz), "%02d" % rep)
            true = read_trees(os.path.join(d, "s_tree.trees"))[0]
            for s in ["0", "25", "50", "100", "250"]:
                for n in ["25", "50", "100", "500"]:
                    for m, pat in files.items():
                        p = os.path.join(d, pat.format(s=s, n=n))
                        if not os.path.exists(p) or os.path.getsize(p) == 0:
                            continue
                        est = read_trees(p)
                        if not est:
                            continue
                        fn, fp, i1, i2 = rf_error(est[0], true)
                        ours = (fn + fp) / (i1 + i2)
                        theirs = pub.get((psiz, float(dl), rep, s, n, m))
                        rows.append([psiz, dl, rep, s, n, m, fn, fp, round(ours, 6), theirs])
                        if theirs is not None:
                            diffs.append(abs(ours - theirs))
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["PSIZ", "DLRT", "REPL", "SQLN", "NGEN", "MTHD", "FN", "FP", "RF_ours", "RF_published"])
    w.writerows(rows)
import statistics
print("trees rescored:", len(rows), "matched to CSV:", len(diffs),
      "exact (|diff|<1e-5):", sum(d < 1e-5 for d in diffs), "max|diff|:", max(diffs))
