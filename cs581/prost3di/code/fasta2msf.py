"""Write a FASTA alignment as GCG MSF (for cross-checking with BAliBASE bali_score)."""
import sys
from bbscore import read_fasta
a = read_fasta(sys.argv[1]); names = list(a); L = len(a[names[0]])
out = ["PileUp", "", "", "", f"   MSF:  {L}  Type: P    Check:  0   ..", ""]
for n in names:
    out.append(f" Name: {n} oo  Len:  {L}  Check:  0  Weight:  10.0")
out += ["", "//", "", ""]
w = max(len(n) for n in names) + 2
for s in range(0, L, 50):
    for n in names:
        seg = a[n][s:s + 50].replace("-", ".")
        out.append(n.ljust(w) + " ".join(seg[i:i + 10] for i in range(0, len(seg), 10)))
    out += ["", ""]
open(sys.argv[2], "w").write("\n".join(out) + "\n")
