"""Mean per-family ASTRAL-Pro quartet support (correct minus each wrong topology) with
standard errors from disjoint blocks, for one 4-taxon pool (heavy-tailed copy numbers
make single big-sample margins unreliable).
usage: python block_se.py POOL.nwk BLOCK OUT.json models"""
import json, os, random, sys
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core

pool, B, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
models = sys.argv[4].split(",")
fams = [l.strip() for l in open(pool) if ";" in l]
spi = {s: i for i, s in enumerate("ABCD")}
TOPO = "(((A,B),C),D);"


def job(arg):
    m, b = arg
    model, par = m.split(":")
    blk = [core.true_rt(f, spi) for f in fams[b * B:(b + 1) * B]]
    if model == "own":
        sc = core.quartet_scores([core.plain(r) for r, _ in blk], [r for r, _ in blk], False, TOPO)
    else:
        rng = random.Random(b * 7919 + len(m))
        err = [core.apply_error(model, float(par), rt, tg, spi, rng) for rt, tg in blk]
        sc = core.quartet_scores([core.encode(r, t) for r, t in err], [r for r, _ in err], True, TOPO)
    d = {t: s for _, t, _, s in sc}
    return m, b, [d["t1"], d["t2"], d["t3"]]


if __name__ == "__main__":
    nb = len(fams) // B
    res = {m: [None] * nb for m in models}
    with mp.Pool(int(os.environ.get("NPROC", "2"))) as p:
        for m, b, v in p.imap_unordered(job, [(m, b) for m in models for b in range(nb)]):
            res[m][b] = v
    summ = {}
    for m, blocks in res.items():
        import statistics as st
        d2 = [(x[0] - x[1]) / B for x in blocks]
        d3 = [(x[0] - x[2]) / B for x in blocks]
        tot = [sum(x) / B for x in blocks]
        se = lambda v: st.stdev(v) / len(v) ** 0.5
        summ[m] = {"per_family_total": st.mean(tot), "c_minus_AD|BC": st.mean(d2), "se2": se(d2),
                   "c_minus_AC|BD": st.mean(d3), "se3": se(d3), "z2": st.mean(d2) / se(d2) if se(d2) else None,
                   "z3": st.mean(d3) / se(d3) if se(d3) else None, "blocks": nb, "block": B}
        print(m, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in summ[m].items()}, flush=True)
    json.dump({"pool": pool, "summary": summ, "blocks": res}, open(out, "w"))
