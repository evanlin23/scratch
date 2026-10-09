"""Sample-complexity curves on the Nute et al. complete ('full') data: 25 ingroup + outgroup,
6 conditions x 20 reps, true and RAxML gene trees, first k genes.
Per species-tree branch: its length in coalescent units (generations / Ne, Ne = 2e5 estimated
by estimate_ne.py) and whether each method recovered it.
Output: /opt/runs/nute_curves.jsonl (restartable).
"""
import sys, os, json, time
from multiprocessing import Pool
import treeswift as ts
sys.path.insert(0, os.path.dirname(__file__))
import stree

ROOT = "/opt/data/nute"
OUT = "/opt/runs/nute_curves.jsonl"
NE = 2.0e5
KS = [10, 25, 50, 100, 250, 500, 1000]


def branch_lengths(true_nwk):
    """bipartition (side without the min label) -> length in CU."""
    t = ts.read_tree_newick(true_nwk)
    allset = frozenset(l.label for l in t.traverse_leaves())
    ref = min(allset)
    desc, out = {}, {}
    for nd in t.traverse_postorder():
        desc[nd] = frozenset([nd.label]) if nd.is_leaf() else frozenset().union(*(desc[c] for c in nd.children))
        s = desc[nd]
        if 1 < len(s) < len(allset) - 1 and nd.parent is not None:
            key = s if ref not in s else allset - s
            out[key] = out.get(key, 0.0) + nd.edge_length / NE
    return out


def run(job):
    h, r, rep, gt = job
    d = f"{ROOT}/25tax-1000gen-0bps-{h}-{r}-full/{rep:02d}"
    true = open(f"{d}/true-species.tre").read().strip()
    bl = branch_lengths(true)
    genes = [stree.strip(g) for g in stree.read_trees(f"{d}/{'true' if gt == 'true' else 'raxml'}-genes.tre")]
    taxa = sorted(stree.bipartitions(true)[1])
    acc = stree.Accumulator(taxa)
    rows, done = [], 0
    for k in KS:
        for g in genes[done:k]:
            acc.add(g)
        done = k
        M = acc.matrix()
        ests, secs = {}, {}
        t0 = time.time(); ests["astrid"] = stree.fastme(M, taxa); secs["astrid"] = time.time() - t0
        ests["njst"] = stree.nj(M, taxa); secs["njst"] = None
        t0 = time.time(); ests["astral"] = stree.astral(genes[:k]); secs["astral"] = time.time() - t0
        for m, e in ests.items():
            eb, _ = stree.bipartitions(e, taxa)
            rec = {",".join(sorted(b)): [round(bl[b], 5), b in eb] for b in bl}
            rows.append(dict(height=h, rate=r, rep=rep, genetrees=gt, k=k, method=m,
                             FN=sum(not v[1] for v in rec.values()) / len(rec), branches=rec, sec=secs[m]))
    return rows


if __name__ == "__main__":
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            x = json.loads(line)
            done.add((x["height"], x["rate"], x["rep"], x["genetrees"]))
    jobs = [(h, r, rep, gt) for rep in range(1, 21) for h in ["10M", "2M", "500K"] for r in ["1E-7", "1E-6"]
            for gt in ["true", "raxml"] if (h, r, rep, gt) not in done]
    print(len(jobs), "jobs", flush=True)
    with Pool(int(os.environ.get("NPROC", 2))) as pool, open(OUT, "a") as fh:
        for rows in pool.imap_unordered(run, jobs):
            for x in rows:
                fh.write(json.dumps(x) + "\n")
            fh.flush()
