"""Follow-up controls on the hard regimes (same seeds as run_exp.py, so the same datasets):
  * saturation handling for the distance baselines: cap factor 2 (main), 1.2, 5, and
    'pcap' (p clipped at 0.75*(1-1/sqrt(k)) before the JC/K2P correction; no inf)
  * FastME-guided controls: Forest components only + GTM(FastME); centroid decomposition + GTM(FastME)
  * threshold-graph-only DCM variant: components of the graph d < m at the selected m
    (this is exactly 'Forest comps only'), reported for FastME guide
usage: python followup.py OUT.jsonl [nproc]
"""
import sys, os, json, tempfile, shutil, warnings, traceback
from multiprocessing import Pool
import numpy as np
import dendropy
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as P  # noqa: E402
from run_exp import key  # noqa: E402
from common import (read_fasta, jc_distance, k2p_distance, cap_distances, tree_splits,  # noqa: E402
                    splits_from_newick, score_tree, _pair_counts)

CELLS = [(100, "U:0.1:0.4", "JC"), (100, "UH:2.0:1.0", "K2P")]
KS = [300, 3000, 100000]
REPS = 10


def pcap_distance(X, model):
    n, k = X.shape
    same, ts = _pair_counts(X)
    eps = 1.0 / np.sqrt(k)
    if model == "JC":
        p = np.minimum(1 - same / k, 0.75 * (1 - eps))
        D = -0.75 * np.log(1 - 4 * p / 3)
    else:
        P_ = ts / k
        Q = np.minimum(1 - same / k - P_, 0.5 * (1 - eps))
        a1 = np.maximum(1 - 2 * P_ - Q, eps * (1 - 2 * Q))
        D = -0.5 * np.log(a1) - 0.25 * np.log(1 - 2 * Q)
    np.fill_diagonal(D, 0)
    return D


def work(args):
    n, k, reg, model, rep = args
    seed = int.from_bytes(f"{n}{k}{reg}{model}".encode(), "little") % 1000003 * 100 + rep
    w = tempfile.mkdtemp(prefix="fu_")
    try:
        rec = dict(n=n, k=k, regime=reg, model=model, rep=rep, key=key(n, k, reg, model, rep))
        tf, aln = P.simulate(n, k, reg, model, seed, w)
        names, X = read_fasta(aln)
        tns = dendropy.TaxonNamespace(names)
        T = tree_splits(dendropy.Tree.get(path=tf, schema="newick", taxon_namespace=tns,
                                          rooting="force-unrooted"), names)
        D0 = jc_distance(X) if model == "JC" else k2p_distance(X)
        Ds = {"cap2": cap_distances(D0, 2.0), "cap1.2": cap_distances(D0, 1.2),
              "cap5": cap_distances(D0, 5.0), "pcap": pcap_distance(X, model)}
        for dn, D in Ds.items():
            for meth, code in (("NJ", "N"), ("FastME", "B")):
                nw = P.fastme(D, names, code, w, meth)
                rec[f"FN_{meth}_{dn}"] = score_tree(splits_from_newick(nw, names), T, n)[0]
        D = Ds["cap2"]
        nj = P.fastme(D, names, "N", w, "nj")
        fme = P.fastme(D, names, "B", w, "fme")
        fo, _, _ = P.best_forest(D0, P.forest_grid(D0, nj, names))
        if fo is None:
            comps = [np.array([i]) for i in range(n)]
            fo = dict(comps=comps, splits=[set() for _ in comps])
        nos = dict(fo, splits=[set() for _ in fo["comps"]])
        for name, f, guide, code in (("FGTM_FastME", fo, fme, "B"), ("CompGTM_FastME", nos, fme, "B")):
            nw, _ = P.forest_gtm(f, guide, names, D, code, "sub", w)
            rec[f"FN_{name}"] = score_tree(splits_from_newick(nw, names), T, n)[0]
        for sz in (25, 50):
            subsets = P.centroid_decomposition(splits_from_newick(fme, names), n, sz)
            dec = dict(comps=subsets, splits=[set() for _ in subsets])
            nw, _ = P.forest_gtm(dec, fme, names, D, "B", "sub", w)
            rec[f"FN_DecGTM_FastME_{sz}"] = score_tree(splits_from_newick(nw, names), T, n)[0]
        # Forest+GTM with the pcap FastME tree as guide and pcap matrix for subsets
        Dp = Ds["pcap"]
        fme_p = P.fastme(Dp, names, "B", w, "fmep")
        nw, _ = P.forest_gtm(fo, fme_p, names, Dp, "B", "sub", w)
        rec["FN_FGTM_FastME_pcap"] = score_tree(splits_from_newick(nw, names), T, n)[0]
        return rec
    except Exception as e:
        return {"key": key(n, k, reg, model, rep), "error": repr(e), "tb": traceback.format_exc()}
    finally:
        shutil.rmtree(w, ignore_errors=True)


def main():
    out = sys.argv[1]
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)["key"] for l in open(out) if l.strip() and "error" not in l}
    tasks = [(n, k, reg, model, rep) for rep in range(REPS) for (n, reg, model) in CELLS for k in KS
             if key(n, k, reg, model, rep) not in done]
    print(len(tasks), "tasks", flush=True)
    with Pool(nproc) as pool, open(out, "a") as fh:
        for r in pool.imap_unordered(work, tasks):
            fh.write(json.dumps(r, default=float) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
