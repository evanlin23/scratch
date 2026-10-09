"""Write the validation report (Markdown) comparing our numbers to the MAGUS paper.

    python -m gcmx.validation_report PUBLISHED.jsonl PAPER_VALUES.json RUNS_DIR PILOT_RESULTS.jsonl > report.md

1. Paper figures vs our FastSP rescoring of the paper's own published alignments.
2. Our MAGUS reruns (paper flags) vs the published MAGUS alignment of the same replicate.
3. Harness check: re-running only the merge on the cached subalignments/backbones
   ("default" pilot variant) vs the full pipeline output it was cached from.
"""

import collections
import glob
import json
import os
import sys

ORDER = ["1000L3", "1000S1", "1000M2", "RNASim", "1000L1", "1000S2", "1000S3", "1000M3", "1000L2", "1000M4",
         "1000M1", "16S.3", "16S.T", "RNASim10K", "16S.B.ALL", "BBA0081", "BBA0067", "BBA0154", "BBA0101",
         "BBA0190", "BBA0134", "BBA0117", "BBA0039", "16S.M"]


def key(dataset, rep):
    return rep.replace("RV100_", "") if dataset == "balibase" else dataset


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def pct(x):
    return "–" if x != x else "{:.1f}".format(100 * x)


def main():
    published_path, paper_path, runs_dir, pilot_path = sys.argv[1:5]
    paper = json.load(open(paper_path))
    rows = [json.loads(l) for l in open(published_path)] if os.path.exists(published_path) else []
    pub = collections.defaultdict(list)
    pub_by_rep = {}
    for r in rows:
        if "avgErr" in r:
            pub[key(r["dataset"], r["rep"]), r["method"]].append(r)
            pub_by_rep[key(r["dataset"], r["rep"]), r["rep"] if r["dataset"] != "balibase" else "R0", r["method"]] = r

    print("## 1. Paper figures vs. our rescoring of the published alignments\n")
    print("Average error (SPFN+SPFP)/2 in %. *Paper* = read off the paper's figures (±0.5); "
          "*rescored* = FastSP on the alignments the authors published (Illinois Data Bank), mean over replicates.\n")
    print("| dataset | reps | MAGUS(Fast) paper | MAGUS(Fast) rescored | PASTA paper | PASTA rescored | MAGUS(Slow) rescored | PASTA(3)+GCM rescored |")
    print("|---|---|---|---|---|---|---|---|")
    for d in ORDER:
        fast = pub.get((d, "MAGUS(Fast)"), [])
        if not fast and d not in paper["error"]:
            continue
        pv = paper["error"].get(d, {})
        print("| {} | {} | {} | {} | {} | {} | {} | {} |".format(
            d, len(fast), pct(pv.get("MAGUS(Fast)", float("nan"))), pct(mean([r["avgErr"] for r in fast])),
            pct(pv.get("PASTA(3)", float("nan"))), pct(mean([r["avgErr"] for r in pub.get((d, "PASTA(3)"), [])])),
            pct(mean([r["avgErr"] for r in pub.get((d, "MAGUS(Slow)"), [])])),
            pct(mean([r["avgErr"] for r in pub.get((d, "PASTA(3)+GCM"), [])]))))

    print("\n## 2. Our MAGUS reruns vs. the published MAGUS alignment of the same replicate\n")
    print("Our runs: current MAGUS (commit 39041fc) with the paper's MAGUS(Fast) flags "
          "(`--maxsubsetsize 0 --maxnumsubsets 25 -r 10 -m 200 --graphbuildhmmextend false`), 4 cores. "
          "Backbone sampling is random, so equal-in-expectation, not identical, results are expected.\n")
    print("| replicate | published MAGUS(Fast) | ours | diff (pts) | published PASTA(3) | our minutes (4 cores) | published minutes |")
    print("|---|---|---|---|---|---|---|")
    diffs = []
    for path in sorted(glob.glob(os.path.join(runs_dir, "*_R*", "prep.json"))):
        ours = json.load(open(path))
        name = os.path.basename(os.path.dirname(path))
        d, rep = name.rsplit("_", 1)
        theirs = pub_by_rep.get((d, rep, "MAGUS(Fast)"))
        pasta = pub_by_rep.get((d, rep, "PASTA(3)"))
        diff = 100 * (ours["avgErr"] - theirs["avgErr"]) if theirs else float("nan")
        if theirs:
            diffs.append(diff)
        print("| {} | {} | {} | {} | {} | {:.0f} | {} |".format(
            name, pct(theirs["avgErr"]) if theirs else "–", pct(ours["avgErr"]),
            "{:+.2f}".format(diff) if theirs else "–", pct(pasta["avgErr"]) if pasta else "–",
            ours["seconds"] / 60, "{:.0f}".format(theirs["published_seconds"] / 60) if theirs and "published_seconds" in theirs else "–"))
    if diffs:
        print("\nMean difference ours − published: {:+.2f} points over {} replicates "
              "(mean |diff| {:.2f}).".format(mean(diffs), len(diffs), mean([abs(x) for x in diffs])))

    print("\n## 3. Harness check: merge-only rerun on cached inputs reproduces the full pipeline\n")
    print("| replicate | full pipeline | merge-only rerun (`default`) | identical error? |")
    print("|---|---|---|---|")
    pilot = {}
    if os.path.exists(pilot_path):
        for r in map(json.loads, open(pilot_path)):
            if r["variant"] == "default" and "avgErr" in r:
                pilot[r["dataset"]] = r
    for path in sorted(glob.glob(os.path.join(runs_dir, "*_R*", "prep.json"))):
        ours = json.load(open(path))
        name = os.path.basename(os.path.dirname(path))
        if name in pilot:
            same = abs(pilot[name]["SPFN"] - ours["SPFN"]) < 1e-9 and abs(pilot[name]["SPFP"] - ours["SPFP"]) < 1e-9
            print("| {} | {} | {} | {} |".format(name, pct(ours["avgErr"]), pct(pilot[name]["avgErr"]), "yes" if same else "NO"))


if __name__ == "__main__":
    main()
