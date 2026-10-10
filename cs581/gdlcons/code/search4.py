"""Random search over 4- and 5-taxon GDL rate configurations for limiting failures of
ASTRAL-Pro under rooting/tagging error models.  Quartet scores come from ASTRAL-Pro3
itself (-C -u 3 on the true species tree), so the 'margin' is
(score of true topology - best wrong score) / total, per internal branch.

usage: python search4.py N_CONFIGS NFAM OUT.jsonl [SEED0]
"""
import json
import multiprocessing as mp
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import (apply_error, encode, plain, quartet_scores, simulate, true_rt)  # noqa: E402
from phylo import parse_newick  # noqa: E402

MODELS = [("true", 0), ("ovl", 0), ("rovl", 0.2), ("rovl", 0.5), ("rtrue", 0.2),
          ("flip", 0.1), ("flip", 0.3), ("d2s", 0.3), ("d2s", 1.0), ("s2d", 0.5), ("own", 0)]

SHAPES = {
    "cat4": "(((A:{a},B:{b})x:{x},C:{c})y:{y},D:{d});",
    "bal4": "((A:{a},B:{b})x:{x},(C:{c},D:{d})y:{y});",
    "cat5": "((((A:{a},B:{b})x:{x},C:{c})y:{y},D:{d})z:{z},E:{e});",
}
LAM = [0, 0, 0.5, 1, 2, 4]
MU = [0, 0.5, 1, 2, 4]
TS = [0.05, 0.2, 0.5, 1.0, 2.0]


def make_config(seed):
    rng = random.Random(seed)
    shape = rng.choice(list(SHAPES))
    names = "abcdexyz"
    lens = {k: rng.choice(TS) for k in names}
    nwk = SHAPES[shape].format(**lens)
    st = parse_newick(nwk)
    rates = {}
    for v in range(len(st.parent)):
        lab = st.label[v] or "root"
        rates[lab] = (rng.choice(LAM), rng.choice(MU))
    root_len = rng.choice([0, 0, 0.5, 1.0])
    return {"seed": seed, "shape": shape, "species": nwk, "rates": rates, "root_len": root_len}


def norm_split(s):
    # '{C}|{D}#{B}|{A}'  -> frozenset of two frozensets (unrooted quartet split)
    a, b = s.split("#")
    side = lambda x: frozenset(y.strip("{}") for y in x.split("|"))
    return frozenset([side(a) | frozenset(), side(b)])


def evaluate(cfg, nfam, seed):
    st = parse_newick(cfg["species"])
    lam = [cfg["rates"][st.label[v] or "root"][0] for v in range(len(st.parent))]
    mu = [cfg["rates"][st.label[v] or "root"][1] for v in range(len(st.parent))]
    try:
        fams, tries, over = simulate(st, lam, mu, nfam, seed, min_species=4,
                                     root_len=cfg["root_len"], cap=3000, max_over=nfam, max_tries=60 * nfam)
    except RuntimeError as e:
        return {"error": str(e)}
    species = sorted(set(st.label[v] for v in st.leaves()))
    spi = {s: i for i, s in enumerate(species)}
    base = [true_rt(f, spi) for f in fams]
    rng = random.Random(seed + 7)
    topo = parse_newick(cfg["species"]).newick()
    res = {}
    for model, par in MODELS:
        if model == "own":
            sc = quartet_scores([plain(r) for r, _ in base], [r for r, _ in base], False, topo)
        else:
            out = [apply_error(model, par, rt, tg, spi, rng) for rt, tg in base]
            sc = quartet_scores([encode(r, t) for r, t in out], [r for r, _ in out], True, topo)
        # rows per branch: t1 is the topology of the given (true) tree
        per = {}
        for node, t, _, score in sc:
            per.setdefault(node, {})[t] = score
        margins = []
        for node, d in per.items():
            tot = sum(d.values()) or 1.0
            margins.append(((d["t1"] - max(d["t2"], d["t3"])) / tot, d["t1"], d["t2"], d["t3"]))
        res["%s:%g" % (model, par)] = margins
    return {"tries": tries, "over": over, "res": res}


def job(args):
    seed, nfam = args
    cfg = make_config(seed)
    r = evaluate(cfg, nfam, seed * 31 + 1)
    r.update(cfg)
    r["nfam"] = nfam
    return r


if __name__ == "__main__":
    n, nfam, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    s0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)["seed"] for l in open(out)}
    todo = [(s, nfam) for s in range(s0, s0 + n) if s not in done]
    with mp.Pool(4) as pool:
        for r in pool.imap_unordered(job, todo):
            with open(out, "a") as f:
                f.write(json.dumps(r) + "\n")
            if "res" in r:
                worst = {k: round(min(m[0] for m in v), 3) for k, v in r["res"].items()}
                print(r["seed"], r["shape"], worst, flush=True)
