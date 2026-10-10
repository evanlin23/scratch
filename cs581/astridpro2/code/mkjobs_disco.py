"""Job list for the DISCO data (100 species, 1000 estimated gene-family trees from 100 bp).
Usage: python mkjobs_disco.py OUT.jsonl JOBS.txt REPS (e.g. 01-10)"""
import json, sys
R = "/opt/data/disco/trees"
CONDS = ["default", "gdl_1e-10_1", "gdl_1e-10_05", "gdl_1e-10_0", "gdl_5e-10_05", "gdl_5e-10_0",
         "gdl_1e-9_1", "gdl_1e-9_05", "gdl_1e-9_0", "ils_1e4", "ils_2e8", "missing_1000"]
a, b = (int(x) for x in sys.argv[3].split("-"))
with open(sys.argv[2], "w") as f:
    for rep in ["%02d" % i for i in range(a, b + 1)]:
        for c in CONDS:
            d = f"{R}/{c}/{rep}"
            key = {"data": "disco", "cond": c, "rep": rep, "sqln": 100, "ngen": 1000}
            f.write(" ".join(["--genes", f"{d}/g_100.trees", "--true", f"{d}/s_tree.trees", "--out", sys.argv[1],
                              "--key", "'" + json.dumps(key) + "'"]) + "\n")
