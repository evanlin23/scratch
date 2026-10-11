"""Why does NeuralNJ do badly on ROSE subsets? Same 50-taxon subsets (true-tree centroid
decomposition), three versions of the data:
  rose      the ROSE true alignment (indels, 1700-4000 columns)
  nogap     the same with every column containing a gap removed (no indel signal, shorter)
  alisim    AliSim GTR+G4 re-simulation on the true subtree, branch lengths and model fitted by
            IQ-TREE on the ROSE alignment (-te true subtree), length 1000, no indels
Usage: python3 diag_nnj.py COND REP OUTDIR NSUB
"""
import json
import os
import subprocess
import sys
import tempfile

import common
from phylo import read_tree, fn_fp
from run_rep_methods import sub_method

cond, rep, outd, nsub = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
d = f"{outd}/{cond}/{rep}"
os.makedirs(d, exist_ok=True)
aln, truef = common.paths(cond, rep)
T = read_tree(truef)
seqs = common.read_fasta(aln)
subs = sorted(common.centroid_decomp(T, 50), key=len, reverse=True)[:nsub]
rows = []
for i, s in enumerate(subs):
    base = f"{d}/s{i}"
    if not os.path.exists(f"{base}.rose.fa"):
        common.write_fasta(seqs, s, f"{base}.rose.fa")
        r = common.read_fasta(f"{base}.rose.fa")
        L = len(r[s[0]])
        keep = [j for j in range(L) if all(r[n][j] != "-" for n in s)]
        with open(f"{base}.nogap.fa", "w") as f:
            for n in s:
                f.write(f">{n}\n{''.join(r[n][j] for j in keep)}\n")
        with open(f"{base}.true.tre", "w") as f:
            f.write(T.copy().restrict(s).to_newick() + "\n")
        with tempfile.TemporaryDirectory() as td:
            # AliSim "mimic" mode: fit GTR+G4 + branch lengths on the fixed true subtree, then simulate
            subprocess.run([common.IQ, "-s", f"{base}.rose.fa", "-te", f"{base}.true.tre", "-m", "GTR+G4", "-T", "1",
                            "--alisim", f"{td}/sim", "--length", "1000", "-af", "fasta", "-seed", str(100 + i),
                            "-pre", f"{td}/fit", "-quiet"], check=True, capture_output=True)
            os.replace(f"{td}/sim.fa", f"{base}.alisim.fa")
    Ts = T.copy().restrict(s)
    for var in ("rose", "nogap", "alisim"):
        a = f"{base}.{var}.fa"
        ncol = len(next(iter(common.read_fasta(a).values())))
        for m in ("FT", "IQ", "BME", "NNJ"):
            o = f"{base}.{var}.{m}.tre"
            if not os.path.exists(o):
                sub_method(m)(a, o)
            fn = fn_fp(read_tree(o), Ts)[0]
            rows.append(dict(cond=cond, rep=rep, subset=i, ntaxa=len(s), variant=var, ncol=ncol, method=m, FN=fn))
            print(cond, rep, i, var, ncol, m, "%.3f" % fn, flush=True)
with open(f"{d}/diag.jsonl", "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
