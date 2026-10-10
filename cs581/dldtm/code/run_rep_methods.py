"""Subset-tree estimators by name (shared by run_rep.py and val_small.py)."""
import os

import common


def sub_method(m):
    if m == "FT":
        return lambda a, o: common.run_fasttree(a, o)
    if m == "IQ":
        return lambda a, o: common.run_iqtree(a, o, threads=1, fast=True)
    if m == "BME":
        return lambda a, o: common.run_fastme_seq(a, o, "B", "J", spr=True)
    if m == "NJ":
        return lambda a, o: common.run_fastme_seq(a, o, "N", "J")
    if m.startswith("PF"):
        import pfdist
        ck = {"PF": f"{common.PF_REPO}/models/pf.ckpt",
              "PFnt": os.environ.get("PFNT_CKPT", "/opt/models/pf_nt.ckpt")}[m]

        def f(a, o):
            names, ss = pfdist.read_fasta(a)
            pfdist.fastme_tree(pfdist.distances(pfdist.load_model(ck), names, ss), names, o)
        return f
    if m == "NNJ":
        import nnj
        return lambda a, o: nnj.tree(a, o)
    raise ValueError(m)
