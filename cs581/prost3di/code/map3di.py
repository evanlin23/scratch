"""Map an alignment of 3Di strings back onto the amino-acid sequences (1:1 residues)."""
import sys
from bbscore import read_fasta
a3 = read_fasta(sys.argv[1]); aa = read_fasta(sys.argv[2])
for n, row in a3.items():
    s = iter(aa[n]); out = "".join("-" if c == "-" else next(s) for c in row)
    assert len(aa[n]) == sum(c != "-" for c in row), n
    print(">" + n); print(out)
