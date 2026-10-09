"""Summarize pilot results.jsonl into a Markdown table.

    python -m gcmx.summarize RESULTS.jsonl [PREP_GLOB] > results.md

Rows: variants; columns: model conditions (replicate dirs are named
CONDITION_Rk). Cells: mean average SP error (SPFN+SPFP)/2 in percent, with
the number of replicates; plus mean alignment-length ratio and runtime.
"""

import collections
import glob
import json
import sys


def condition(dataset):
    return dataset.rsplit("_R", 1)[0]


def main():
    rows = [json.loads(line) for line in open(sys.argv[1])]
    if len(sys.argv) > 2:
        for path in glob.glob(sys.argv[2]):
            prep = json.load(open(path))
            name = prep["outdir"].rstrip("/").rsplit("/", 1)[1]
            rows.append(dict(prep, dataset=name, variant="MAGUS full pipeline (prep)"))
    rows = [r for r in rows if "avgErr" in r]

    conds = sorted({condition(r["dataset"]) for r in rows})
    variants = list(dict.fromkeys(r["variant"] for r in rows))
    cell = collections.defaultdict(list)
    for r in rows:
        cell[r["variant"], condition(r["dataset"])].append(r)

    print("| variant | " + " | ".join(conds) + " | len ratio | sec |")
    print("|---|" + "---|" * (len(conds) + 2))
    for v in variants:
        parts, ratios, secs = [], [], []
        for c in conds:
            rs = cell[v, c]
            if rs:
                parts.append("{:.2f} (n={})".format(100 * sum(r["avgErr"] for r in rs) / len(rs), len(rs)))
                ratios += [r["LenEst"] / r["LenRef"] for r in rs if r.get("LenRef")]
                secs += [r["seconds"] for r in rs]
            else:
                parts.append("–")
        ratio = "{:.2f}".format(sum(ratios) / len(ratios)) if ratios else "–"
        sec = "{:.0f}".format(sum(secs) / len(secs)) if secs else "–"
        print("| {} | {} | {} | {} |".format(v, " | ".join(parts), ratio, sec))

    # paired comparison against default, per replicate
    base = {r["dataset"]: r["avgErr"] for r in rows if r["variant"] == "default"}
    print("\nPaired difference vs default (percentage points, negative = better), wins/ties/losses:\n")
    print("| variant | mean diff | W/T/L |")
    print("|---|---|---|")
    for v in variants:
        diffs = [100 * (r["avgErr"] - base[r["dataset"]]) for r in rows
                 if r["variant"] == v and r["dataset"] in base and v != "default"]
        if diffs:
            w = sum(d < -0.05 for d in diffs)
            l = sum(d > 0.05 for d in diffs)
            print("| {} | {:+.2f} | {}/{}/{} |".format(v, sum(diffs) / len(diffs), w, len(diffs) - w - l, l))


if __name__ == "__main__":
    main()
