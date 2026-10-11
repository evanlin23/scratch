"""Job lists. Usage: python mkjobs.py {train|heldout} OUT.jsonl JOBS.txt METHODS [THREADS] [fmrfs|disco|both]"""
import json, sys
F = "/opt/data/fmrfs"
DS = "/opt/data/disco"
DLS = ("0.0000000001", "0.0000000002", "0.0000000005")
PSS = ("10000000", "50000000")


def fmrfs(reps):
    for rep in reps:
        for dl in DLS:
            for ps in PSS:
                for sq in (25, 100):
                    d = f"{F}/ntaxa-100.dlrate-{dl}.psize-{ps}/{rep}"
                    yield ({"data": "fmrfs", "cond": f"dl{dl[-1]}e-10_ps{ps[0]}e7", "rep": rep, "sqln": sq},
                           f"{d}/ft-sqlen-{sq}-ngen-100.trees", f"{d}/s_tree.trees")


def disco(pairs):
    for c, rep in pairs:
        yield ({"data": "disco", "cond": c, "rep": rep, "sqln": 100},
               f"{DS}/sequences/{c}/{rep}/ft_100.trees", f"{DS}/trees/{c}/{rep}/s_tree.trees")


def sets(which):
    if which == "train":
        return list(fmrfs(["01", "02", "03"])), list(disco([(c, r) for r in ("01", "02") for c in ("default", "gdl_1e-10_05")]))
    return (list(fmrfs(["%02d" % i for i in range(4, 11)])),
            list(disco([(c, r) for r in ("03", "04") for c in ("default", "gdl_1e-10_05")] +
                       [(c, r) for c in ("gdl_5e-10_05", "ils_2e8") for r in ("01", "02", "03", "04")])))


if __name__ == "__main__":
    which, out, jobs, methods = sys.argv[1:5]
    th = sys.argv[5] if len(sys.argv) > 5 else "1"
    part = sys.argv[6] if len(sys.argv) > 6 else "both"
    fm, dc = sets(which)
    items = (fm if part in ("fmrfs", "both") else []) + (dc if part in ("disco", "both") else [])
    with open(jobs, "w") as f:
        for key, g, t in items:
            f.write(json.dumps({"genes": g, "true": t, "key": key, "out": out, "methods": methods, "threads": th}) + "\n")
