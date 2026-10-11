"""Build 'true' (experimental-structure) 3Di strings for BAliBASE sequences.

For each sequence, pick the chain of its PDB entry (first 4 characters of the name) whose structure
sequence aligns best to the BAliBASE sequence (local alignment, BLOSUM62), and transfer the 3Di state of
each aligned residue; residues without a structure counterpart get 'X' (score 0 in the 3Di matrix).
Usage: python3 true3di.py STRUCT_AA.fa STRUCT_3DI.fa OUTDIR set1.aa.fa [set2.aa.fa ...]
Prints per-set coverage (fraction of residues with a transferred state, and identity of matched pairs).
"""
import os, sys, collections
from Bio import Align
from Bio.Align import substitution_matrices
from bbscore import read_fasta

def read_raw(path):
    d, n = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            n = line[1:].split()[0]; d[n] = []
        else:
            d[n].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}

saa, s3 = read_raw(sys.argv[1]), read_raw(sys.argv[2])
byid = collections.defaultdict(list)
for k in saa:
    byid[k[:4].lower()].append(k)
al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=-10, extend_gap_score=-1)
out = sys.argv[3]; os.makedirs(out, exist_ok=True)
for f in sys.argv[4:]:
    seqs = read_fasta(f); res = {}; cov = []; ids = []
    for n, s in seqs.items():
        best = None
        for c in byid.get(n[:4].lower(), []):
            a = al.align(s, saa[c])[0]
            if best is None or a.score > best[0]:
                best = (a.score, c, a)
        t = ["X"] * len(s)
        if best:
            _, c, a = best; same = tot = 0
            for (qs, qe), (ts, te) in zip(*a.aligned):
                for i, j in zip(range(qs, qe), range(ts, te)):
                    t[i] = s3[c][j]; tot += 1; same += s[i] == saa[c][j]
            ids.append(same / max(tot, 1))
        cov.append(sum(x != "X" for x in t) / len(s))
        res[n] = "".join(t)
    b = os.path.basename(f).replace(".aa.fa", "")
    with open(os.path.join(out, b + ".3di.fa"), "w") as fh:
        for n, t in res.items():
            fh.write(f">{n}\n{t}\n")
    print(b, f"{sum(cov)/len(cov):.3f}", f"{min(cov):.3f}", f"{sum(ids)/max(len(ids),1):.3f}")
