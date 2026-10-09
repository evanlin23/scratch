"""DEV-only screen of correction variants (odd replicates) on the hardest missing-data settings."""
import sys, os, json
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
import stree, md_exp

def run(job):
    h, r, pat, ng, rep, gt = job
    genes, true = md_exp.load(h, r, pat, ng, rep, gt)
    taxa = sorted(stree.bipartitions(true)[1]); ph = stree.estimate_p(genes, len(taxa))
    out = []
    for mode in ["plain", "ht1", "plugin"]:
        M, C = stree.distance_matrix(genes, taxa, mode=mode, p=ph)
        for tm, fn in [("fastme", stree.fastme), ("nj", stree.nj)]:
            out.append(dict(height=h, rate=r, pattern=pat, ngenes=ng, rep=rep, genetrees=gt, mode=mode, tree=tm,
                            FN=stree.fn_fp(true, fn(M, taxa))[0]))
    return out

if __name__ == "__main__":
    jobs = [(h, r, pat, ng, rep, gt) for rep in range(1, 21, 2) for h in ["10M", "2M", "500K"] for r in ["1E-7", "1E-6"]
            for gt in ["true", "raxml"] for pat, ng in [("rand-60", 200), ("rand-60", 1000), ("sim-70", 1000)]]
    with Pool(2) as p, open("/opt/runs/dev_variants.jsonl", "w") as fh:
        for rows in p.imap_unordered(run, jobs):
            for x in rows: fh.write(json.dumps(x) + "\n")
