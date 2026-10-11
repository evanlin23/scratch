"""Positive control: NeuralNJ vs FastTree / IQ-TREE / BME / NJ on NeuralNJ's own released test
set (GTR+I+G AliSim, 50 taxa, 256/512/1024 sites; data_gen/data/test in the NeuralNJ repo).
Checks that our CPU harness reproduces the paper's in-distribution accuracy.
Usage: python3 indomain_nnj.py OUTDIR N_PER_LENGTH OUT_JSONL
"""
import glob
import json
import os
import subprocess
import sys

import common
from phylo import read_tree, fn_fp
from run_rep_methods import sub_method

outd, nper, outj = sys.argv[1], int(sys.argv[2]), sys.argv[3]
os.makedirs(outd, exist_ok=True)
rows = []
for L in (256, 512, 1024):
    phys = sorted(p for p in glob.glob(f"/opt/src/NeuralNJ/data_gen/data/test/len{L}/taxa50/*.phy")
                  if os.path.exists(p[:-4] + ".tre"))[:nper]
    for phy in phys:
        b = f"{outd}/{os.path.basename(phy)[:-4]}"
        T = read_tree(phy[:-4] + ".tre")
        lines = [l for l in open(phy).read().split("\n")[1:] if l.strip()]
        with open(b + ".fa", "w") as f:
            for l in lines:
                p = l.split()
                f.write(f">{p[0]}\n{''.join(p[1:])}\n")
        for m in ("FT", "IQ", "BME", "NJ", "NNJ"):
            o = f"{b}.{m}.tre"
            if not os.path.exists(o):
                if m == "NNJ":   # feed the original PHYLIP (keeps NeuralNJ's own taxon naming)
                    subprocess.run([common.os.environ.get("NNJ_PY", "/opt/mm/root/envs/nnj/bin/python"),
                                    os.path.join(os.path.dirname(os.path.abspath(__file__)), "nnj_driver.py"), phy, o],
                                   check=True, capture_output=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
                else:
                    sub_method(m)(b + ".fa", o)
            rows.append(dict(kind="indomain", length=L, file=os.path.basename(phy), method=m,
                             FN=fn_fp(read_tree(o), T)[0]))
            print(L, os.path.basename(phy), m, "%.3f" % rows[-1]["FN"], flush=True)
with open(outj, "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
