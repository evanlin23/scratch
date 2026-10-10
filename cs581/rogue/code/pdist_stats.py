import sys, random, statistics as st
sys.path.insert(0, "/home/user/scratch/cs581/code")
from gcmx import fasta
a = fasta.read(sys.argv[1] + "/true.fasta"); R = set(open(sys.argv[1] + "/rogues.txt").read().split())
C = [n for n in a if n not in R]
def p(x, y):
    s = t = 0
    for c1, c2 in zip(a[x], a[y]):
        if c1 != "-" and c2 != "-":
            t += 1; s += c1 != c2
    return s / max(t, 1)
rng = random.Random(0)
for lab, grp in (("rogue", sorted(R)[:15]), ("clean", rng.sample(C, 15))):
    nn = [min(p(x, y) for y in rng.sample(C, 300) if y != x) for x in grp]
    print(lab, "nearest-of-300 p: median %.3f" % st.median(nn), "seqlen", st.median(len(a[x].replace('-', '')) for x in grp))
