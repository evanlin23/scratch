"""Run the identifiability-aware post-processor (idclass.py) on every stored DecoDiPhy run of the
authors' E1 simulations (DecoDiPhy-Data/biotrees/*/*_k*_s*_noise*/all_rounds.json).
For each run: flags of the TRUE configuration and of DecoDiPhy's FINAL fit, the equivalence class of
the fit on its edges, and accuracy of DecoDiPhy's point vs the canonical class representative.
usage: python analyze_data.py DATA_DIR OUT.jsonl [nproc]"""
import sys, os, glob, json, re
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from multiprocessing import Pool
from idclass import URTree


def ekey(pl):
    return tuple(sorted(pl[:2]))


def run(f):
    try:
        R = URTree(open(f + "/pruned_tree.trees").read().strip().split("\n")[0])
        rounds = json.load(open(f + "/all_rounds.json"))
        tr = rounds[0]
        fins = [r for r in rounds if r["rounds"] == "final"]
        fin = fins[-1] if fins else rounds[-1]
        out = dict(run=f.split("/")[-1], tree=f.split("/")[-2])
        mm = re.search(r"_k(\d+)_s(\d+)_(noise\d)(?:_l(\d+))?", out["run"])
        out.update(k=int(mm.group(1)), seed=int(mm.group(2)), noise=mm.group(3), l=mm.group(4) or "-")
        res = {}
        for name, r in (("true", tr), ("fit", fin)):
            pls = [R.placement(a, x) for a, x in zip(r["anchors"], r["x"])]
            fl = R.flags(pls, r["x"])
            fl.pop("W")
            res[name] = (r, pls, fl)
            out[name + "_k"] = len(pls)
            for kk, v in fl.items():
                out[f"{name}_{kk}"] = bool(v)
        rt, plt, _ = res["true"]
        rf, plf, _ = res["fit"]
        Et = {ekey(p) for p in plt}; Ef = {ekey(p) for p in plf}
        out["jaccard"] = len(Et & Ef) / len(Et | Ef)
        out["same_edges"] = Et == Ef
        truth = [(pl, pq) for pl, pq in zip(plt, rt["p"])]
        fitm = [(pl, pq) for pl, pq in zip(plf, rf["p"])]
        out["emd_fit"] = R.emd(truth, fitm)
        C = R.eq_class(plf, rf["p"], rf["y"])
        if C is None:
            out["class_ok"] = False
            return out
        out["class_ok"] = True
        width = C["hi"] - C["lo"]
        out["p_width_max"] = float(width.max())
        out["p_width_sum"] = float(width.sum())
        out["nonid_class"] = bool(width.max() > 1e-6)
        canon = [(pl, pq) for pl, pq in zip(C["canon_pls"], C["canon_p"])]
        out["emd_canon"] = R.emd(truth, canon)
        ve = [R.emd(truth, list(zip(vp, vq))) for vq, vp in C["vertices"]]
        out["emd_class_min_vertex"] = float(min(ve)); out["emd_class_max_vertex"] = float(max(ve))
        # can some placement slide to a node within the class (edge ambiguity)?
        out["slide_to_node"] = bool(np.any(C["lo"] < -1)) if False else None
        if out["same_edges"]:
            # abundance L1 error, matched by edge
            pt = {ekey(pl): pq for pl, pq in zip(plt, rt["p"])}
            ordr = [ekey(pl) for pl in plf]
            ptv = np.array([pt[e] for e in ordr])
            out["pL1_fit"] = float(np.abs(np.array(rf["p"]) - ptv).sum())
            out["pL1_canon"] = float(np.abs(C["canon_p"] - ptv).sum())
            out["truth_in_p_intervals"] = bool(np.all((ptv >= C["lo"] - 1e-3) & (ptv <= C["hi"] + 1e-3)))
            out["fit_p_equal_truth"] = bool(np.all(np.abs(np.array(rf["p"]) - ptv) <= 1e-3))
            # best abundance error achievable inside the class (interval clipping lower bound)
            out["pL1_class_lb"] = float(np.sum(np.maximum(0, np.maximum(C["lo"] - ptv, ptv - C["hi"]))))
        return out
    except Exception as e:
        return dict(run=f, error=repr(e))


if __name__ == "__main__":
    runs = sorted(glob.glob(sys.argv[1] + "/*/*_k*_s*_noise*"))
    npr = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    with Pool(npr) as P, open(sys.argv[2], "w") as fo:
        for i, o in enumerate(P.imap_unordered(run, runs, chunksize=4)):
            fo.write(json.dumps(o) + "\n")
            if i % 200 == 0:
                print(i, len(runs), flush=True)
