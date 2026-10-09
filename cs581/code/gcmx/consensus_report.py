"""Significance analysis of the consensus-MSA test (pre-registered before results).

    python -m gcmx.consensus_report RESULTS.jsonl [RESULTS.jsonl ...] > REPORT.md

Primary hypotheses, per scenario, paired by replicate, two-sided Wilcoxon
signed-rank test on the average error (SPFN+SPFP)/2:
  H1  consensus vs the "central" input (the input the method picks as primary;
      what a user could choose without a reference alignment)
  H2  consensus vs the best input of each replicate (an oracle; conservative)
Secondary: consensus vs each fixed input, and per-condition breakdowns.
Differences are in percentage points (negative = consensus better);
W/T/L uses a tie band of |diff| <= 0.05 points.
"""

import collections
import json
import sys

from scipy.stats import wilcoxon

LABELS = {"gcm.txt": "MAGUS(Fast)", "gcm_slow.txt": "MAGUS(Slow)", "pasta_align.txt": "PASTA(3)",
          "pasta_3_gcm_align.txt": "PASTA(3)+GCM", "pasta_1_align.txt": "PASTA(1)", "pasta_4_align.txt": "PASTA(4)"}


def condition(rep):
    ds = rep.split("/")[0]
    return "BAliBASE" if ds == "balibase" else ds


def stats(diffs):
    n = len(diffs)
    w = sum(d < -0.05 for d in diffs)
    l = sum(d > 0.05 for d in diffs)
    try:
        p = wilcoxon(diffs).pvalue if n >= 6 and any(diffs) else float("nan")
    except ValueError:
        p = float("nan")
    return n, sum(diffs) / n, "{}/{}/{}".format(w, n - w - l, l), p


def fmt_p(p):
    return "–" if p != p else ("<1e-4" if p < 1e-4 else "{:.3g}".format(p))


def main():
    rows = {}
    for path in sys.argv[1:]:
        for line in open(path):
            r = json.loads(line)
            if "consensus" in r and "avgErr" in r["consensus"]:
                rows[r["replicate"], r["scenario"]] = r  # last write wins (reruns)
    by_scen = collections.defaultdict(list)
    for (rep, scen), r in sorted(rows.items()):
        by_scen[scen].append(r)

    print("# Consensus-MSA significance test\n")
    print("Inputs are the MAGUS paper's published alignments of each replicate; the consensus is "
          "`gcmx.consensus` (fixed design: most-central input as primary, ~13-sequence similarity groups, "
          "all inputs as GCM evidence). Error = (SPFN+SPFP)/2 in %.\n")
    for scen, rs in by_scen.items():
        inputs = rs[0]["inputs"]
        cons = [100 * r["consensus"]["avgErr"] for r in rs]
        best = [100 * min(r["input_err"].values()) for r in rs]
        central = [100 * r["input_err"][r["central_input"]] for r in rs]
        print("## Scenario: {} ({} replicates)\n".format(scen, len(rs)))
        print("| alignment | mean error | consensus − it (pts) | W/T/L | Wilcoxon p |")
        print("|---|---|---|---|---|")
        print("| **consensus** | {:.2f} | | | |".format(sum(cons) / len(cons)))
        for label, base in (("H1: central input", central), ("H2: best input (oracle)", best)):
            n, mean, wtl, p = stats([c - b for c, b in zip(cons, base)])
            print("| {} | {:.2f} | {:+.2f} | {} | {} |".format(label, sum(base) / n, mean, wtl, fmt_p(p)))
        for f in inputs:
            base = [100 * r["input_err"][f] for r in rs]
            n, mean, wtl, p = stats([c - b for c, b in zip(cons, base)])
            print("| {} | {:.2f} | {:+.2f} | {} | {} |".format(LABELS.get(f, f), sum(base) / n, mean, wtl, fmt_p(p)))
        picked = collections.Counter(LABELS.get(r["central_input"], r["central_input"]) for r in rs)
        print("\nCentral input chosen: " + ", ".join("{} {}×".format(k, v) for k, v in picked.most_common()) + "\n")

        print("| condition | n | consensus | central | best input | cons − central | W/T/L | p | cons − best | W/T/L | p |")
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        groups = collections.defaultdict(list)
        for r, c, b, ce in zip(rs, cons, best, central):
            groups[condition(r["replicate"])].append((c, b, ce))
        for g in sorted(groups):
            cs, bs, ces = zip(*groups[g])
            n1, m1, w1, p1 = stats([c - x for c, x in zip(cs, ces)])
            n2, m2, w2, p2 = stats([c - x for c, x in zip(cs, bs)])
            print("| {} | {} | {:.2f} | {:.2f} | {:.2f} | {:+.2f} | {} | {} | {:+.2f} | {} | {} |".format(
                g, n1, sum(cs) / n1, sum(ces) / n1, sum(bs) / n1, m1, w1, fmt_p(p1), m2, w2, fmt_p(p2)))
        print()


if __name__ == "__main__":
    main()
