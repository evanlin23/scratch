"""Re-estimate FastMulRFS-data gene trees WITH branch support, from the published alignments.

For each condition / replicate / sequence length L in {25, 100}: the first 100 genes listed in
genes-with-gt3species.txt (the same genes, in the same order, as the published RAxML trees
g_trees-raxml-sqlen-L.trees, which were estimated on the first L sites), truncated to the first L
sites, FastTree 2 -nt -gtr -gamma (SH-like local support in [0,1] on internal edges).
Output: <rep dir>/ft-sqlen-L-ngen-100.trees.  Usage: python estimate_fmrfs.py [NPROC]
"""
import glob, os, subprocess, sys, tempfile
from multiprocessing import Pool

R = "/opt/data/fmrfs"
FT = "/opt/mm/root/envs/gdl/bin/FastTree"
NG = 100


def read_phy(p):
    with open(p) as f:
        f.readline()
        for line in f:
            x = line.split()
            if len(x) == 2:
                yield x[0], x[1]


def one(job):
    d, L = job
    out = f"{d}/ft-sqlen-{L}-ngen-{NG}.trees"
    if os.path.exists(out):
        return out
    genes = [l.strip() for l in open(f"{d}/genes-with-gt3species.txt") if l.strip()][:NG]
    trees = []
    with tempfile.TemporaryDirectory(dir="/opt/tmp") as td:
        for g in genes:
            fa = os.path.join(td, "a.fa")
            with open(fa, "w") as f:
                for n, s in read_phy(f"{d}/{int(g):04d}.phy"):
                    f.write(f">{n}\n{s[:L]}\n")
            r = subprocess.run([FT, "-nt", "-gtr", "-gamma", "-quiet", "-nopr", fa], capture_output=True, text=True, check=True)
            trees.append(r.stdout.strip())
    with open(out + ".tmp", "w") as f:
        f.write("\n".join(trees) + "\n")
    os.rename(out + ".tmp", out)
    return out


if __name__ == "__main__":
    jobs = [(d, L) for d in sorted(glob.glob(f"{R}/ntaxa-100.*/[0-9][0-9]")) for L in (25, 100)]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as p:
        for i, o in enumerate(p.imap_unordered(one, jobs)):
            print(i, o, flush=True)
