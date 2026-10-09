"""Limiting quartet sums with seed-to-seed variation for one configuration.
Usage: python margin.py TREE RATES NFAM SEEDS... (prints JSON per seed)"""
import itertools
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe import run  # noqa: E402

tree, rates_s, nfam = sys.argv[1], sys.argv[2], int(sys.argv[3])
rates = {k: tuple(v) for k, v in json.loads(rates_s).items()}
leaves = sorted(set(c for c in tree if c.isupper()))
qs = list(itertools.combinations(leaves, 4))
for seed in sys.argv[4:]:
    r = run(tree, rates, (0, 0), nfam, int(seed), 2, astral=False, quartets=qs,
            modes=[("multi", "mean", False), ("pro", "mean", True), ("pro", "mean", False)])
    r["seed"] = seed
    print(json.dumps(r), flush=True)
