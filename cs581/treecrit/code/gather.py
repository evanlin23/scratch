# AI-assisted (Claude), exploration code for CS581 project
"""Collect tree rows (FastTree nRF vs true tree) and FastSP rows for every bank replicate in reps.tsv.

    python3 gather.py   -> data/trees.jsonl, data/aln_logged.jsonl  (one row per key x method)

Rows are read with `git show origin/<branch>:cs581/<results_dir>/{trees,aln,vote_aln,vote}.jsonl`.
Tree method names are normalised to: true, magus, es4, recipe, es3, hardf (gg linsi&fftns2-op3),
vote_hard, vote_hard-bb, vote_soft, vote_soft-bb, vote_soft2, vote_soft4, vote_hard_mask (masked output),
split_magus, split_recipe.
"""
import json
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, "..", "data")

TREE_NAME = {"hard": "hardf", "vote_hard+mask": "vote_hard_mask", "vote_hard_mask_masked": "vote_hard_mask"}
GG_NAME = {"linsi": "magus", "wsoft0.03:linsi&fftns2#es4": "recipe", "linsi#es3": "es3",
           "linsi&fftns2-op3": "hardf", "linsi#es4": "es4"}


def show(br, path):
    r = subprocess.run(["git", "-C", REPO, "show", "origin/{}:{}".format(br, path)], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else ""


def rows(text):
    for l in text.splitlines():
        l = l.strip()
        if l:
            try:
                yield json.loads(l)
            except ValueError:
                pass


def main():
    trees, alns = [], []
    for line in open(os.path.join(HERE, "reps.tsv")):
        if line.startswith("#") or not line.strip():
            continue
        key, br, bank, res, src = line.rstrip("\n").split("\t")
        name = key.rsplit("_", 1)[0]
        base = "cs581/" + res
        seen = {}
        for r in rows(show(br, base + "/trees.jsonl")):
            if (r.get("dataset") or r.get("rep")) != name:
                continue
            m = TREE_NAME.get(r["method"], r["method"])
            if m == "magus_shuf":
                continue
            seen[m] = {"key": key, "dataset": name, "method": m, "RF": r["RF"], "FN": r["FN"], "FP": r["FP"]}
        trees += seen.values()
        aseen = {}
        for f in ("aln.jsonl", "vote_aln.jsonl", "vote.jsonl"):
            for r in rows(show(br, base + "/" + f)):
                if r.get("rep") != name or "avgErr" not in r:
                    continue
                v = r["variant"]
                if "B" in r and r["B"] not in (10, None):
                    continue
                m = GG_NAME.get(v) if "B" not in r else ("magus" if v == "magus" else "es4" if v == "es4" else
                                                          "vote_" + v.replace("+mask", "_mask"))
                if m is None:
                    continue
                if m in aseen and not ("B" in r and m == "es4"):  # trees used vote.py's es4 (vote/es4/out.fasta)
                    continue
                aseen[m] = ({"key": key, "method": m, "variant": v, **{k: r[k] for k in
                                     ("SPFN", "SPFP", "avgErr", "TC", "LenEst", "LenRef")},
                                     **({"masked_cols": r["masked_cols"]} if "masked_cols" in r else {})})
        alns += aseen.values()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "trees.jsonl"), "w") as f:
        for r in trees:
            f.write(json.dumps(r) + "\n")
    with open(os.path.join(OUT, "aln_logged.jsonl"), "w") as f:
        for r in alns:
            f.write(json.dumps(r) + "\n")
    print(len(trees), "tree rows,", len(alns), "logged alignment rows")
    import collections
    c = collections.Counter(r["method"] for r in trees)
    print(dict(c))
    have = {(r["key"], r["method"]) for r in alns}
    miss = [(r["key"], r["method"]) for r in trees if r["method"] not in ("true", "split_magus", "split_recipe")
            and (r["key"], r["method"]) not in have]
    print("tree rows without a logged SP row:", miss)


if __name__ == "__main__":
    main()
