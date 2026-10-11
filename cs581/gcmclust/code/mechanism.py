"""es4 gain per clustering method: mean over reps of err(es4:X:p) - err(raw:X:p), proteins vs DNA/RNA."""
import sys
import numpy as np
import summarize as S

rows = S.load()
which = sys.argv[1] if len(sys.argv) > 1 else "train"
reps = [r for r in S.SPLIT if (which == "all" or S.SPLIT[r] == which) and r in rows]
pairs = sorted({v.split(":", 1)[1] for r in reps for v in rows[r]})
print("| method:param | es4 − raw, proteins (n) | es4 − raw, DNA/RNA (n) | raw vs MAGUS, proteins | raw vs MAGUS, DNA/RNA |")
print("|---|---|---|---|---|")
for mp in pairs:
    out = []
    for sel in (S.PROT, lambda r: not S.PROT(r)):
        d = [rows[r]["es4:" + mp]["err"] - rows[r]["raw:" + mp]["err"] for r in reps if sel(r)
             and "err" in rows[r].get("es4:" + mp, {}) and "err" in rows[r].get("raw:" + mp, {})]
        out.append("{:+.2f} ({})".format(np.mean(d), len(d)) if d else "–")
    for sel in (S.PROT, lambda r: not S.PROT(r)):
        d = [rows[r]["raw:" + mp]["err"] - rows[r][S.BASE]["err"] for r in reps if sel(r)
             and "err" in rows[r].get("raw:" + mp, {})]
        out.append("{:+.2f}".format(np.mean(d)) if d else "–")
    print("| {} | {} |".format(mp, " | ".join(out)))
