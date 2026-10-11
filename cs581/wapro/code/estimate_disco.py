"""Re-estimate DISCO-data gene-family trees WITH support: FastTree 2 -nt -gtr -gamma (SH-like local support)
on the published 100 bp alignments of all 1000 families. Families with < 3 sequences get the trivial tree.
Output: /opt/data/disco/sequences/<cond>/<rep>/ft_100.trees.  Usage: python estimate_disco.py NPROC COND,COND REPS"""
import os, subprocess, sys, tempfile
from multiprocessing import Pool

R = "/opt/data/disco/sequences"
FT = "/opt/mm/root/envs/gdl/bin/FastTree"


def gene(p):
    seqs = []
    with open(p) as f:
        f.readline()
        for line in f:
            x = line.split()
            if len(x) == 2:
                seqs.append(x)
    if len(seqs) < 3:
        return "(" + ",".join(n for n, _ in seqs) + ");"
    with tempfile.NamedTemporaryFile("w", suffix=".fa", dir="/opt/tmp", delete=False) as f:
        for n, s in seqs:
            f.write(f">{n}\n{s}\n")
    try:
        r = subprocess.run([FT, "-nt", "-gtr", "-gamma", "-quiet", "-nopr", f.name], capture_output=True, text=True, check=True)
    finally:
        os.unlink(f.name)
    return r.stdout.strip()


def one(d):
    out = f"{d}/ft_100.trees"
    if os.path.exists(out):
        return out
    files = sorted(x for x in os.listdir(f"{d}/100") if x.endswith(".phy"))
    trees = [gene(f"{d}/100/{x}") for x in files]
    with open(out + ".tmp", "w") as f:
        f.write("\n".join(trees) + "\n")
    os.rename(out + ".tmp", out)
    return out


if __name__ == "__main__":
    a, b = (int(x) for x in sys.argv[3].split("-"))
    jobs = [f"{R}/{c}/{r:02d}" for r in range(a, b + 1) for c in sys.argv[2].split(",")]
    with Pool(int(sys.argv[1])) as p:
        for i, o in enumerate(p.imap_unordered(one, jobs)):
            print(i, o, flush=True)
