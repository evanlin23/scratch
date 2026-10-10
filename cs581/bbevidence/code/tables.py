"""Curated report tables from the result files.   python3 tables.py RESULTS_DIR > results/tables.md"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize import load, per_rep, table  # noqa: E402

d = sys.argv[1]
res = load(d, "results")
PROT = {r for r in res if r.startswith("BBA")}
NUC = {r for r in res if not r.startswith("BBA")}

GROUPS = [
    ("Single-tool backbones (10, MAGUS's sequence sets)", ["clustalo", "fftns2", "famsa", "clustalo-iter"]),
    ("Amount and diversity of evidence", ["linsi~5", "clustalo~5", "clustalo@s1", "clustalo@n20", "clustalo@n30",
                                          "linsi+clustalo", "linsi~5+clustalo^5"]),
    ("Consensus: pair intersections (same sequence sets)", ["linsi&clustalo", "linsi&fftns2", "linsi&fftns2-op3",
                                                            "linsi&famsa", "clustalo&fftns2", "linsi&clustalo&fftns2",
                                                            "linsi&clustalo+clustalo@s1&fftns2@s1",
                                                            "clustalo&fftns2+clustalo@s1&fftns2@s1+clustalo@s2&fftns2@s2"]),
    ("Consensus: column masks on MAGUS's own L-INS-i backbones (MAFFT only, no new alignments)",
     ["linsi|cons0.3", "linsi|cons0.5", "linsi|cons0.7", "linsi|cons0.8", "linsi|agree-clustalo-0.5"]),
    ("Union plus consistency mask", ["linsi+clustalo|cons0.5", "linsi+clustalo|cons0.7", "linsi+linsi&clustalo",
                                     "linsi+linsi&fftns2", "linsi+linsi&clustalo&fftns2"]),
    ("MAFFT settings and gappy-column masks", ["linsi-op3@h", "fftns2@h", "linsi@h&fftns2@h", "ginsi@h", "ginsi-ul8@h",
                                               "ginsi-ul4@h", "linsi|gap0.5", "linsi|gap0.7", "linsi|gap0.9",
                                               "clustalo|gap0.7"]),
]
for label, reps in (("Protein (BAliBASE)", PROT), ("Nucleotide (RNASim, 16S.M, ROSE)", NUC)):
    print("## " + label)
    for title, vs in GROUPS:
        print(table(res, vs, reps, title))
    print()
    print("#### Per replicate: control error (%), then Δ vs control (points)")
    key = ["clustalo", "linsi+clustalo", "linsi&clustalo", "linsi&fftns2", "linsi&fftns2-op3", "clustalo&fftns2",
           "linsi|cons0.5", "linsi|cons0.7", "linsi|cons0.8", "linsi+clustalo|cons0.5", "linsi+clustalo|cons0.7"]
    print(per_rep(res, key, reps))
    print()
