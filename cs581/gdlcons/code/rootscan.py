"""How often does ASTRAL-Pro3's own min-duplication rooting break a 4-taxon caterpillar
whose correct-root overlap tagging is consistent?  Random configurations from the
prevalence scan with feasible copy numbers; single-sample margins at NFAM families.
usage: python rootscan.py N NFAM OUT.jsonl [NPROC]"""
import json, os, random, sys
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core
from phylo import parse_newick

def feasible(b):
    return (b['y'][0] - b['y'][1]) * b['y'][2] <= 3.5 and all((b[k][0] - b[k][1]) * b[k][2] <= 3 for k in 'ABCx') \
        and max(b[k][0] * b[k][2] for k in 'ABCxy') <= 16

def job(arg):
    i, r, nfam = arg
    b = r['br']
    nwk = "(((A:%g,B:%g)x:%g,C:%g)y:%g,D:%g);" % (b['A'][2], b['B'][2], b['x'][2], b['C'][2], b['y'][2],
                                                b['D'][2] if 'D' in b else 1)
    st = parse_newick(nwk)
    lam = [0.0] * len(st.parent); mu = [0.0] * len(st.parent)
    for v in range(len(st.parent)):
        if st.label[v] in b:
            lam[v], mu[v] = b[st.label[v]][0], b[st.label[v]][1]
    try:
        fams, tries, over = core.simulate(st, lam, mu, nfam, 1000 + i, min_species=4, cap=4000, max_over=nfam, max_tries=400 * nfam)
    except RuntimeError as e:
        return {"i": i, "br": b, "error": str(e)}
    spi = {s: j for j, s in enumerate("ABCD")}
    base = [core.true_rt(f, spi) for f in fams]
    out = {"i": i, "br": b, "pred_ovl_margin": r['ovl_margin'], "nfam": nfam}
    for m in ("own", "ovl"):
        if m == "own":
            sc = core.quartet_scores([core.plain(x) for x, _ in base], [x for x, _ in base], False, st.newick())
        else:
            err = [core.apply_error("ovl", 0, rt, tg, spi, random.Random(i)) for rt, tg in base]
            sc = core.quartet_scores([core.encode(x, t) for x, t in err], [x for x, _ in err], True, st.newick())
        d = {core.split_of(lab): s for _, _, lab, s in sc}
        c, w1, w2 = d.get("AB|CD", 0), d.get("AC|BD", 0), d.get("AD|BC", 0)
        out[m] = (c - max(w1, w2)) / max(c + w1 + w2, 1e-9)
    return out

if __name__ == "__main__":
    n, nfam, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    rows = [json.loads(l) for l in open('/opt/runs/gdlcons/prevalence.jsonl')]
    cand = [r for r in rows if r['ovl_margin'] is not None and r['O'] > 1e-3 and r['ovl_margin'] > 0 and feasible(r['br'])]
    rng = random.Random(5)
    sel = rng.sample(cand, n)
    with mp.Pool(nproc) as p:
        for r in p.imap_unordered(job, [(i, r, nfam) for i, r in enumerate(sel)]):
            with open(out, "a") as f:
                f.write(json.dumps(r) + "\n")
            print(r.get("i"), r.get("pred_ovl_margin"), r.get("ovl"), r.get("own"), r.get("error", ""), flush=True)
