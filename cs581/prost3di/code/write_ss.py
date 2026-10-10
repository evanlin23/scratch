"""Overwrite the 3Di track (db_ss) of a Foldseek database with strings from a FASTA file.

Usage: python3 write_ss.py DBPREFIX new3di.fa
Sequences are matched by name through DBPREFIX.lookup; lengths must equal the AA lengths.
"""
import sys
from bbscore import read_fasta
db, fa = sys.argv[1], sys.argv[2]
names = {}
for line in open(db + ".lookup"):
    k, n = line.split("\t")[:2]
    names[int(k)] = n
aa_len = {}
for line in open(db + ".index"):
    k, off, ln = map(int, line.split("\t"))
    aa_len[k] = ln - 2
new = read_fasta(fa)
data = bytearray(); idx = []
for k in sorted(names):
    s = new[names[k]].upper()
    assert len(s) == aa_len[k], (names[k], len(s), aa_len[k])
    rec = s.encode() + b"\n\0"
    idx.append(f"{k}\t{len(data)}\t{len(rec)}\n"); data += rec
open(db + "_ss", "wb").write(bytes(data))
open(db + "_ss.index", "w").write("".join(idx))
