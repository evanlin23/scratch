"""Cheap alignment statistics not reported by FastSP (cached).

    python alnstats.py OUT.jsonl ALN [ALN ...]

Per alignment: number of columns, gap fraction, and indel events = maximal internal
gap runs summed over sequences (terminal gaps excluded).
"""
import json
import os
import re
import sys

GAP = re.compile(r"[-.]+")


def read(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def stats(path):
    s = read(path)
    L = len(next(iter(s.values())))
    runs = gaps = 0
    for x in s.values():
        core = x.strip("-.")
        runs += len(GAP.findall(core))
        gaps += x.count("-") + x.count(".")
    return {"ncol": L, "gapfrac": gaps / (L * len(s)), "indels": runs}


if __name__ == "__main__":
    out = sys.argv[1]
    done = {json.loads(l)["aln"] for l in open(out)} if os.path.exists(out) else set()
    for p in sys.argv[2:]:
        if p in done or not os.path.exists(p):
            continue
        r = {"aln": p}
        r.update(stats(p))
        with open(out, "a") as f:
            f.write(json.dumps(r) + "\n")
