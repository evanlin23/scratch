"""Trimmed queue (see round 3 note below) (written after a throughput check, before looking at accuracy):
 1. FastMulRFS data, 100 genes, reps 01-10: all methods.
 2. DISCO data, 1000 genes, reps 01-05: all methods except DupLoss-2 (wQFM-GDL only reps 01-03).
 3. FastMulRFS data, 500 genes, reps 01-05: all methods except wQFM-GDL and DupLoss-2.
Each line: bench.py args incl. --methods."""
import json, sys
out_fm, out_di, jobs = sys.argv[1], sys.argv[2], sys.argv[3]
FAST = "astrid-multi,astrid-pro,astrid-pro-r0,astrid-pro-s,astrid-disco,asteroid"
SLOW = "astral-pro3,disco-astral,fastmulrfs"
L = []
def fm(ngen, rep, methods):
    for dl in ("0.0000000001", "0.0000000002", "0.0000000005"):
        for ps in ("10000000", "50000000"):
            for sq in (25, 100):
                d = f"/opt/data/fmrfs/ntaxa-100.dlrate-{dl}.psize-{ps}/{rep}"
                key = {"data": "fmrfs", "cond": f"dl{dl[-1]}e-10_ps{ps[0]}e7", "rep": rep, "sqln": sq, "ngen": ngen}
                L.append(f"--genes {d}/g_trees-raxml-sqlen-{sq}.trees --true {d}/s_tree.trees --ngen {ngen} --out {out_fm} --key '{json.dumps(key)}' --methods {methods}")
CONDS = ["default", "gdl_1e-10_1", "gdl_1e-10_05", "gdl_1e-10_0", "gdl_5e-10_05", "gdl_5e-10_0",
         "gdl_1e-9_1", "gdl_1e-9_05", "gdl_1e-9_0", "ils_1e4", "ils_2e8", "missing_1000"]
def di(rep, methods):
    for c in CONDS:
        d = f"/opt/data/disco/trees/{c}/{rep}"
        key = {"data": "disco", "cond": c, "rep": rep, "sqln": 100, "ngen": 1000}
        L.append(f"--genes {d}/g_100.trees --true {d}/s_tree.trees --out {out_di} --key '{json.dumps(key)}' --methods {methods}")
# Round 3 (throughput): DupLoss-2 dropped after reps 01-04 (53/1/0 losses); wQFM-GDL kept on
# FastMulRFS reps 01-04 and DISCO reps 01-02 only.
for i in range(1, 6):
    di("%02d" % i, FAST + "," + SLOW + (",wqfm-gdl" if i <= 2 else ""))
for i in range(5, 11):
    fm(100, "%02d" % i, FAST + "," + SLOW)
for i in range(1, 6):
    fm(500, "%02d" % i, FAST + ",astral-pro3")
open(jobs, "w").write("\n".join(L) + "\n")
print(len(L))
