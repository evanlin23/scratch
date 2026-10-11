"""Generate one simulated replicate: Yule-like model tree with short internal branches,
AliSim GTR+G alignment (true alignment, no indels).
Usage: python3 simgen.py OUTDIR N INTERNAL_MEAN REP [SITES]
"""
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../gtm/code"))
from sim import model_tree  # noqa: E402

IQ = "/opt/mm/root/envs/bio/bin/iqtree3"
out, n, imean, rep = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
sites = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
os.makedirs(out, exist_ok=True)
if not os.path.exists(f"{out}/aln.fa"):
    rng = random.Random(7919 * rep + n)
    open(f"{out}/true.tre", "w").write(model_tree("yule", n, rng, imean, 0.1) + "\n")
    subprocess.run([IQ, "--alisim", f"{out}/aln", "-t", f"{out}/true.tre", "-m",
                    "GTR{1.0,3.0,1.0,1.0,3.0}+F{0.3,0.2,0.2,0.3}+G4{0.5}", "--length", str(sites), "-af", "fasta",
                    "-seed", str(rep + 17), "-redo"], check=True, capture_output=True)
    for f in os.listdir(out):
        if f.startswith("aln") and f != "aln.fa":
            os.remove(f"{out}/{f}")
