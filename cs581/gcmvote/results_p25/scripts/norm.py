# AI-assisted (Claude), exploration code for CS581 project
import sys,hashlib
def rd(p):
    d={};n=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): n=l[1:].split()[0]; d[n]=[]
        elif n: d[n].append(l)
    return {k:''.join(v) for k,v in d.items()}
for p in sys.argv[1:]:
    d=rd(p); print(hashlib.md5(''.join(k+d[k] for k in sorted(d)).encode()).hexdigest()[:10], len(d), len(next(iter(d.values()))), p)
