import math, sys
from famR import margin
for lamT in [0.5, 1, 2, 3, 4, 6, math.inf]:
    row = []
    for a in [0.5, 0.2, 0.1, 0.05, 0.01]:
        bestc = None
        for c in [0.9, 0.7, 0.5, 0.3, 0.2, 0.1]:
            if c <= a: continue
            m, r = margin(a, a, c, lamT)
            tot = r["O"] + r["AB"] + r["AC"] + r["BC"]
            if bestc is None or m / tot < bestc[0]:
                bestc = (m / tot, c)
        row.append("a=b=%g:%.3f(c=%g)" % (a, bestc[0], bestc[1]))
    print("lamT=%s" % lamT, " ".join(row))
