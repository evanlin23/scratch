"""Summaries + paired tests. usage: python3 analyze.py <tag> [<tag> ...] > out.md
Condition = case id without the replicate suffix. Paired comparisons are by case."""
import sys, json, os, collections
from scipy.stats import wilcoxon

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
TIE = 0.25  # percentage points of RF: |diff| <= TIE counts as a tie

PAIRS = [  # (new, baseline)
    ("dc:astral3:100:astral3", "astral3"),
    ("dc:mrlft:100:astral3", "astral3"),
    ("dc:mrlft:200:astral3", "astral3"),
    ("dc:mrlft:100:astral3", "mrlft"),
    ("dc:mrlft:200:astral3", "mrlft"),
    ("mrlft", "astral3"),
    ("astral4", "astral3"),
    ("tqmc", "astral3"),
    ("dc:scs:100:astral3", "scs"),
    ("dc:scs:500:astral3", "scs"),
    ("dc:scs:100:astral3", "astral3"),
    ("dc:scs:500:astral3", "astral3"),
    ("scs", "astral3"),
    ("dc:astral3:100:astral3", "scs"),
    ("dc:scs:200:astral3", "scs"),
    ("dc:true:100:astral3", "astral3"),
    ("dc:dc:mrlft:200:astral3:200:astral3", "dc:mrlft:200:astral3"),
]


def cond(case):
    return case.rsplit("_r", 1)[0]


def load(tags):
    rows = {}
    for tag in tags:
        for line in open(os.path.join(RES, f"{tag}.jsonl")):
            r = json.loads(line)
            if "rf" in r:
                rows[(cond(r["case"]), r["case"], r["method"])] = r  # last one wins
    return rows


def main():
    rows = load(sys.argv[1:])
    conds = sorted(set(k[0] for k in rows), key=lambda c: (len(c), c))
    meths = sorted(set(k[2] for k in rows), key=lambda m: (m.startswith("dc"), m))
    print("### Mean RF error (%) / mean FN (%) / mean wall (s) / mean peak RSS (MB), n replicates\n")
    print("| condition | method | n | RF % | FN % | FP % | wall s | peak MB |")
    print("|---|---|---|---|---|---|---|---|")
    for c in conds:
        for m in meths:
            v = [r for (cc, _, mm), r in rows.items() if cc == c and mm == m]
            if not v:
                continue
            f = lambda k: sum(x[k] or 0 for x in v) / len(v)
            print(f"| {c} | {m} | {len(v)} | {100*f('rf'):.2f} | {100*f('fn'):.2f} | {100*f('fp'):.2f} | "
                  f"{f('wall'):.1f} | {f('peak_mb'):.0f} |")
    print(f"\n### Paired comparisons (RF %, new - baseline; negative = new is better). "
          f"W/T/L = new wins/ties/losses, tie band |diff| <= {TIE} pp; two-sided Wilcoxon signed-rank\n")
    print("| condition | new | baseline | n | mean dRF (pp) | W/T/L (RF) | p (RF) | mean dFN (pp) | W/T/L (FN) | p (FN) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for c in conds + ["ALL"]:
        for new, base in PAIRS:
            cols = []
            for key in ("rf", "fn"):
                d = []
                for (cc, case, mm), r in rows.items():
                    if mm != new or (c != "ALL" and cc != c):
                        continue
                    b = rows.get((cc, case, base))
                    if b:
                        d.append(100 * (r[key] - b[key]))
                if len(d) < 3:
                    break
                w = sum(x < -TIE for x in d); l = sum(x > TIE for x in d); t = len(d) - w - l
                try:
                    p = wilcoxon(d).pvalue if any(x != 0 for x in d) else 1.0
                except ValueError:
                    p = float("nan")
                cols.append((len(d), f"{sum(d)/len(d):+.2f} | {w}/{t}/{l} | {p:.2g}"))
            if len(cols) < 2:
                continue
            print(f"| {c} | {new} | {base} | {cols[0][0]} | {cols[0][1]} | {cols[1][1]} |")


if __name__ == "__main__":
    main()
