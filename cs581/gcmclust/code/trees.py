"""FastTree (-lg -gamma) on merged alignments of the AliSim replicates; nRF vs the true tree (dendropy).

    python3 trees.py OUT.jsonl REP:VARIANT [...]     (VARIANT 'true' = the true alignment)
"""
import json
import os
import subprocess
import sys

import dendropy
from dendropy.calculate import treecompare

W = "/opt/work/gcmclust/reps"
out = sys.argv[1]
done = {(r["rep"], r["variant"]) for r in map(json.loads, open(out))} if os.path.exists(out) else set()
for job in sys.argv[2:]:
    rep, var = job.split(":", 1)
    if (rep, var) in done:
        continue
    aln = os.path.join(W, rep, "true.fasta") if var == "true" else os.path.join(
        W, rep, "variants", var.replace(":", "_"), "out.fasta")
    if var == "raw:mcl:4" and not os.path.exists(aln):
        aln = os.path.join(W, rep, "ctrl", "out.fasta")
    tre = aln + ".ft.nwk"
    with open(tre, "w") as f:
        subprocess.run(["FastTree", "-lg", "-gamma", "-quiet", "-nopr", aln], stdout=f, stderr=subprocess.DEVNULL,
                       check=True)
    k, r = rep.split("_R")
    tns = dendropy.TaxonNamespace()
    t_true = dendropy.Tree.get(path="/opt/data/sim/{}/R{}/tree.nwk".format(k, r), schema="newick", taxon_namespace=tns)
    t_est = dendropy.Tree.get(path=tre, schema="newick", taxon_namespace=tns)
    for t in (t_true, t_est):
        t.is_rooted = False
        t.encode_bipartitions()
    rf = treecompare.symmetric_difference(t_true, t_est)
    n = len(tns)
    row = {"rep": rep, "variant": var, "nRF": round(rf / (2 * (n - 3)), 4)}
    with open(out, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
