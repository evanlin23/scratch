"""FAMSA and MAFFT --auto on the MAGUS paper's ROSE 1000M2 and RNASim-1000 replicates (1 thread),
to compare a fast aligner's rank across simulators with the published MAGUS/PASTA errors.

    python real_sims.py OUT.jsonl [reps]
"""
import json, os, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))
from gcmx import score, fasta  # noqa: E402
SETS = {"1000M2": "/opt/data/Datasets/ROSE/1000M2/R{}/rose.aln.true.fasta",
        "RNASim1000": "/opt/data/Datasets/RNASim/1000/R{}/true_align.txt"}
out, reps = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 5
done = {(r["set"], r["rep"], r["method"]) for r in map(json.loads, open(out))} if os.path.exists(out) else set()
for rep in range(reps):
    for s, pat in SETS.items():
        ref = pat.format(rep)
        with tempfile.TemporaryDirectory() as tmp:
            un = os.path.join(tmp, "u.fa")
            seqs = fasta.ungap(fasta.upper(fasta.read(ref)))
            fasta.write({k: v.replace("U", "T") for k, v in seqs.items()}, un)  # RNASim is RNA; FAMSA needs DNA letters
            for m, cmd in (("famsa", ["/opt/mm/root/envs/bio/bin/famsa", "-t", "1", un, os.path.join(tmp, "famsa.fa")]),
                           ("mafft-auto", ["mafft", "--auto", "--thread", "1", un])):
                if (s, rep, m) in done:
                    continue
                o = os.path.join(tmp, m + ".fa")
                with open(o if m == "mafft-auto" else os.devnull, "w") as f:
                    subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, check=True)
                refT = os.path.join(tmp, "ref.fa")
                fasta.write({k: v.upper().replace("U", "T") for k, v in fasta.read(ref).items()}, refT)
                sc = score.fastsp(refT, o)
                with open(out, "a") as f:
                    f.write(json.dumps({"set": s, "rep": rep, "method": m, "avgErr": sc["avgErr"], "SPFN": sc["SPFN"], "SPFP": sc["SPFP"]}) + "\n")
