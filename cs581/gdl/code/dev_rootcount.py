"""Exploratory (dev replicates 01-03 only): ASTRID-Pro with vs without counting the
inferred gene-tree root. Usage: python dev_rootcount.py out.jsonl"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import methods as M  # noqa: E402
from phylo import read_trees, rf_error  # noqa: E402

base = "/opt/data/gdl/fmrfs"
with open(sys.argv[1], "w") as f:
    for dl in ["0.0000000001", "0.0000000002", "0.0000000005"]:
        for ps in ["10000000", "50000000"]:
            for rep in ["01", "02", "03"]:
                d = os.path.join(base, "ntaxa-100.dlrate-%s.psize-%s" % (dl, ps), rep)
                true = read_trees(os.path.join(d, "s_tree.trees"))[0]
                sp = sorted(true.label[v] for v in true.leaves())
                idx = {s: i for i, s in enumerate(sp)}
                for sq in ["25", "100"]:
                    allg = read_trees(os.path.join(d, "g_trees-raxml-sqlen-%s.trees" % sq))
                    for n in [25, 100, 500]:
                        tg = [M.root_and_tag(g, M.simphy_species, idx) for g in allg[:n]]
                        out = {"dl": dl, "ps": ps, "rep": rep, "sqln": sq, "ngen": n}
                        for name, cr in [("pro", False), ("pro-countroot", True)]:
                            per = [M.gene_distances(rt, t, ls, len(sp), mode="pro", count_root=cr) for rt, t, ls in tg]
                            D, _ = M.average_matrix(per)
                            out[name] = rf_error(M.fastme_tree(D, sp), true)[0]
                        f.write(json.dumps(out) + "\n"); f.flush()
