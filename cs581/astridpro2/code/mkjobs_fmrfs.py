"""Job list for the FastMulRFS data: 3 DL rates x 2 psize x 10 reps x sqlen {25,100} x ngen {100,500}."""
import json, os, sys
R = "/opt/data/fmrfs"
out = sys.argv[1]
H = os.path.dirname(os.path.abspath(__file__))
with open(sys.argv[2], "w") as f:
    for ngen in (100, 500):
        for rep in ["%02d" % i for i in range(1, 11)]:
            for dl in ("0.0000000001", "0.0000000002", "0.0000000005"):
                for ps in ("10000000", "50000000"):
                    for sq in (25, 100):
                        d = f"{R}/ntaxa-100.dlrate-{dl}.psize-{ps}/{rep}"
                        key = {"data": "fmrfs", "cond": f"dl{dl[-1]}e-10_ps{ps[0]}e7", "rep": rep, "sqln": sq, "ngen": ngen}
                        f.write(" ".join(["--genes", f"{d}/g_trees-raxml-sqlen-{sq}.trees", "--true", f"{d}/s_tree.trees",
                                          "--ngen", str(ngen), "--out", out, "--key", "'" + json.dumps(key) + "'"]) + "\n")
