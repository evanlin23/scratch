"""Model trees + AliSim alignments for the simulator-dependence pilot.

    python simulate.py OUTDIR --shape bd|balanced|caterpillar --n 500 --sigma 0.3
        --height H --indel 0.05 --size geo:MEAN|pow:A/MAX --seed S

Tree: ultrametric topology+node heights (birth-death, perfectly balanced, or
caterpillar), then each branch length multiplied by an independent log-normal
rate exp(N(-sigma^2/2, sigma^2)) (mean 1; uncorrelated lineage-specific rate
heterogeneity = clock violation), then rescaled so the mean root-to-tip
distance equals H (expected substitutions per site). sigma = 0 is a strict clock.

Sequences: AliSim (IQ-TREE 3), HKY (kappa 2, base freqs .3/.2/.2/.3) + G4 (alpha 1),
root length 1000, insertion rate = deletion rate = INDEL/2 relative to the
substitution rate, indel lengths geometric (mean MEAN) or power-law (Zipf a,
truncated at MAX). Writes OUTDIR/{model.tre,true.fasta,unaligned.fasta,params.json}.
"""
import argparse
import json
import math
import os
import random
import subprocess

IQTREE = "/opt/mm/root/envs/bio/bin/iqtree3"


class Node:
    __slots__ = ("children", "height", "length", "name")

    def __init__(self, height, children=(), name=None):
        self.children, self.height, self.name, self.length = list(children), height, name, 0.0


def bd_tree(n, rng, birth=1.0, death=0.5):
    """Reconstructed birth-death tree on n tips: node heights from the
    coalescent point process (Gernhard 2008 / Stadler), uniformly random ranked topology."""
    # CPP: n-1 i.i.d. node depths with the BD conditioned-on-n density, inverse-CDF sampling
    # (Gernhard 2008, with lambda=birth, mu=death, conditioned on origin ~ uniform prior).
    lam, mu = birth, death
    r = lam - mu
    depths = []
    for _ in range(n - 1):
        u = rng.random()
        # inverse CDF of H with F(t) = (1 - e^{-rt}) / (1 - (mu/lam) e^{-rt}) * ... use simple form:
        # for the "uniform origin" prior the node depths satisfy t = (1/r) ln((lam - mu u)/(lam (1-u)))
        depths.append(math.log((lam - mu * u) / (lam * (1 - u))) / r)
    # coalescent point process: tips in a line, node between tip i and i+1 at depth depths[i]
    return cpp_build(depths)


def cpp_build(depths):
    """Build the tree whose i-th internal node (between leaves i and i+1) has height depths[i]."""
    n = len(depths) + 1
    leaves = [Node(0.0, name="t{}".format(i)) for i in range(n)]

    def build(lo, hi):  # leaves lo..hi inclusive
        if lo == hi:
            return leaves[lo]
        k = max(range(lo, hi), key=lambda i: depths[i])
        return Node(depths[k], [build(lo, k), build(k + 1, hi)])
    import sys
    sys.setrecursionlimit(10000)
    return build(0, n - 1)


def balanced_tree(n, H):
    leaves = iter(range(n))

    def build(m, h):
        if m == 1:
            return Node(0.0, name="t{}".format(next(leaves)))
        return Node(h, [build(m // 2, h - step), build(m - m // 2, h - step)])
    depth = math.ceil(math.log2(n))
    step = H / depth
    return build(n, H)


def caterpillar_tree(n, rng):
    # node heights from the same BD node-height distribution as bd_tree, sorted: the caterpillar
    # attaches one tip per internal node
    tmp = bd_tree(n, rng)
    hs = sorted(internal_heights(tmp))
    node = Node(0.0, name="t0")
    for i, h in enumerate(hs):
        node = Node(h, [node, Node(0.0, name="t{}".format(i + 1))])
    return node


def internal_heights(t):
    out, stack = [], [t]
    while stack:
        x = stack.pop()
        if x.children:
            out.append(x.height)
            stack.extend(x.children)
    return out


def assign_lengths(t, sigma, H, rng):
    stack = [t]
    while stack:
        x = stack.pop()
        for c in x.children:
            mult = math.exp(rng.gauss(-sigma * sigma / 2, sigma)) if sigma > 0 else 1.0
            c.length = (x.height - c.height) * mult
            stack.append(c)
    # rescale to mean root-to-tip H
    rtt, stack = [], [(t, 0.0)]
    while stack:
        x, d = stack.pop()
        if not x.children:
            rtt.append(d)
        for c in x.children:
            stack.append((c, d + c.length))
    s = H / (sum(rtt) / len(rtt))
    stack = [t]
    while stack:
        x = stack.pop()
        x.length *= s
        stack.extend(x.children)
    rtt = [r * s for r in rtt]
    m = sum(rtt) / len(rtt)
    return {"rtt_cv": (sum((r - m) ** 2 for r in rtt) / len(rtt)) ** 0.5 / m, "rtt_max_over_min": max(rtt) / min(rtt)}


def newick(t):
    def rec(x):
        if not x.children:
            return "{}:{:.6g}".format(x.name, x.length)
        return "({}):{:.6g}".format(",".join(rec(c) for c in x.children), x.length)
    return "(" + ",".join(rec(c) for c in t.children) + ");"


def size_spec(s):
    kind, val = s.split(":")
    if kind == "geo":
        return "GEO{{{:.6g}}}".format(float(val))  # AliSim GEO{mean}
    if kind == "pow":
        return "POW{{{}}}".format(val)
    raise ValueError(s)


def read_fasta(p):
    seqs, name = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--shape", default="bd")
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--sigma", type=float, default=0.3)
    ap.add_argument("--height", type=float, default=1.0)
    ap.add_argument("--indel", type=float, default=0.05)
    ap.add_argument("--size", default="geo:4")
    ap.add_argument("--length", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--model", default="HKY{2}+F{0.3/0.2/0.2/0.3}")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    if a.shape == "bd":
        t = bd_tree(a.n, rng)
    elif a.shape == "balanced":
        t = balanced_tree(a.n, 1.0)
    elif a.shape == "caterpillar":
        t = caterpillar_tree(a.n, rng)
    else:
        raise ValueError(a.shape)
    tstats = assign_lengths(t, a.sigma, a.height, rng)
    tre = os.path.join(a.out, "model.tre")
    with open(tre, "w") as f:
        f.write(newick(t) + "\n")
    prefix = os.path.join(a.out, "alisim")
    sz = size_spec(a.size)
    cmd = [IQTREE, "--alisim", prefix, "-t", tre, "-m", a.model, "--length", str(a.length),
           "--indel", "{0:.6g},{0:.6g}".format(a.indel / 2), "--indel-size", "{0},{0}".format(sz),
           "--seed", str(a.seed), "-af", "fasta", "--no-unaligned", "-redo"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    aln = read_fasta(prefix + ".fa")
    aln = {k: v.upper() for k, v in aln.items() if k.startswith("t")}
    # drop all-gap columns
    names = list(aln)
    keep = [i for i in range(len(aln[names[0]])) if any(aln[k][i] != "-" for k in names)]
    with open(os.path.join(a.out, "true.fasta"), "w") as f, open(os.path.join(a.out, "unaligned.fasta"), "w") as g:
        for k in names:
            s = "".join(aln[k][i] for i in keep)
            f.write(">{}\n{}\n".format(k, s))
            g.write(">{}\n{}\n".format(k, s.replace("-", "")))
    for ext in (".fa", ".log"):
        if os.path.exists(prefix + ext):
            os.remove(prefix + ext)
    json.dump({**vars(a), **tstats, "alisim_cmd": " ".join(cmd)}, open(os.path.join(a.out, "params.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
