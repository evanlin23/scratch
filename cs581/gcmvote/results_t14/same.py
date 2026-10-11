# AI-assisted (Claude), exploration code for CS581 project
"""Compare two FASTA alignments as {name: row} (order-insensitive)."""
import sys
def rd(p):
    d={};k=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): k=l[1:].strip(); d[k]=[]
        elif k: d[k].append(l)
    return {k:''.join(v) for k,v in d.items()}
a=rd(sys.argv[1]); b=rd(sys.argv[2]); print(a==b, len(next(iter(a.values()))), len(next(iter(b.values()))))
