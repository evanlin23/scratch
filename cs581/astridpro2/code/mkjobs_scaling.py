"""Job list for runtime scaling: (a) #genes on DISCO gtrees_10000_l1 (100 species, estimated trees from
100 bp), reps 01-02, ngen 100..10000; (b) #taxa on induced subsets of species_1000 rep 01 (mk_scaling.py)."""
import json, sys
out, jobs = sys.argv[1], sys.argv[2]
R = "/opt/data/disco/trees/gtrees_10000_l1"
with open(jobs, "w") as f:
    for k in (50, 100, 200, 500, 1000):
        d = f"/opt/data/scaling/taxa_{k}"
        key = {"data": "scale_taxa", "cond": f"taxa_{k}", "rep": "01", "ngen": 500}
        f.write(f"--genes {d}/genes.trees --true {d}/s_tree.trees --out {out} --key '{json.dumps(key)}'\n")
    for rep in ("01", "02"):
        for n in (100, 300, 1000, 3000, 10000):
            key = {"data": "scale_genes", "cond": "gtrees_10000_l1", "rep": rep, "ngen": n}
            f.write(f"--genes {R}/{rep}/g_100.trees --true {R}/{rep}/s_tree.trees --ngen {n} --out {out} --key '{json.dumps(key)}'\n")
