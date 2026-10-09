"""Random search over per-branch GDL rates on 4- and 5-taxon species trees for cases where a
method's limiting distance matrix gives the wrong tree (FastME FN>0 with many families).
Usage: python search.py SEED NCONFIG NFAM out.jsonl"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe import run  # noqa: E402

seed, ncfg, nfam, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
rng = random.Random(seed)
SHAPES = {
    "cat4": ("(((A:{A},B:{B})x:{x},C:{C})y:{y},D:{D})r;", "ABCDxy"),
    "bal4": ("((A:{A},B:{B})x:{x},(C:{C},D:{D})y:{y})r;", "ABCDxy"),
    "cat5": ("((((A:{A},B:{B})x:{x},C:{C})y:{y},D:{D})z:{z},E:{E})r;", "ABCDExyz"),
}
LAM = [0, 0.3, 1, 2, 4]
MU = [0, 0.3, 1, 2, 4]
with open(out, "a") as f:
    for i in range(ncfg):
        shape = rng.choice(list(SHAPES))
        tmpl, nodes = SHAPES[shape]
        bl = {n: round(rng.choice([0.1, 0.3, 1.0, 2.0]), 2) for n in nodes}
        rates = {n: (rng.choice(LAM), rng.choice(MU)) for n in nodes}
        tree = tmpl.format(**bl)
        try:
            r = run(tree, rates, (0, 0), nfam, seed * 100000 + i, 2, astral=False)
        except Exception as e:  # noqa
            continue
        r.update({"shape": shape, "tree": tree, "rates": rates})
        f.write(json.dumps(r) + "\n")
        f.flush()
