"""Job lists. Usage: python mkjobs.py EXP OUT.jsonl JOBS.txt
EXP fmrfs : FastMulRFS 6 conds x reps 01-10, 100 genes, levels true (+tags) / est100 / est25
EXP disco : DISCO 100-species conditions, reps given by env REPS (default 01-03), levels true (+tags) / est100"""
import json, os, sys
exp, out, jobs = sys.argv[1:4]
MAIN = "astrid-multi,astrid-pro,astrid-disco,disco-astral,astral-pro3,asteroid,wqfm-gdl"
TT = "astrid-pro-tt,astrid-disco-tt,disco-astral-tt"
HYB = "apro3-fast,apro3-guide-fast,apro3-guide"
L = []
def job(genes, true, key, methods, tagged=None, ngen=0, threads=1):
    a = ["--genes", genes, "--true", true, "--out", out, "--key", "'" + json.dumps(key) + "'", "--methods", methods,
         "--threads", str(threads)]
    if tagged: a += ["--tagged", tagged]
    if ngen: a += ["--ngen", str(ngen)]
    L.append(" ".join(a))
if exp == "fmrfs":
    for rep in ["%02d" % i for i in range(1, 11)]:
        for dl in ("0.0000000001", "0.0000000002", "0.0000000005"):
            for ps in ("10000000", "50000000"):
                c = f"dl{dl[-1]}e-10_ps{ps[0]}e7"
                d = f"/opt/data/fmrfs/ntaxa-100.dlrate-{dl}.psize-{ps}/{rep}"
                p = f"/opt/data/prep/fmrfs/{c}/{rep}"
                k = {"data": "fmrfs", "cond": c, "rep": rep, "ngen": 100}
                job(p + "/true.trees", d + "/s_tree.trees", dict(k, level="true"), MAIN + "," + TT, tagged=p + "/tt.trees")
                for sq in (100, 25):
                    job(f"{d}/g_trees-raxml-sqlen-{sq}.trees", d + "/s_tree.trees", dict(k, level="est%d" % sq),
                        MAIN + "," + HYB, ngen=100)
elif exp == "disco":
    a, b = (int(x) for x in os.environ.get("REPS", "01-03").split("-"))
    conds = os.environ.get("CONDS", "default,gdl_1e-9_1,gdl_5e-10_0,gdl_1e-10_0,ils_2e8,ils_1e4").split(",")
    meth = os.environ.get("METHODS", MAIN)
    for rep in ["%02d" % i for i in range(a, b + 1)]:
        for c in conds:
            d = f"/opt/data/disco/trees/{c}/{rep}"
            p = f"/opt/data/prep/disco/{c}/{rep}"
            k = {"data": "disco", "cond": c, "rep": rep, "ngen": 1000}
            for lev in os.environ.get("LEVELS", "true,est100").split(","):
                if lev == "true":
                    job(d + "/g_true.trees", d + "/s_tree.trees", dict(k, level="true"), meth + "," + TT, tagged=p + "/tt.trees")
                else:
                    job(f"{d}/g_{lev[3:]}.trees", d + "/s_tree.trees", dict(k, level=lev), meth)
open(jobs, "w").write("\n".join(L) + "\n")
print(len(L), "jobs")
