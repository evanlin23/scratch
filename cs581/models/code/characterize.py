"""Summary statistics of a true alignment (+ optional model tree), used to match AliSim to ROSE.

    python characterize.py TRUE.fasta [TREE.newick]  -> JSON

Stats: n taxa, mean ungapped length, true-alignment length, gap fraction,
mean pairwise p-distance (sampled pairs, over columns where both have letters),
gap-run length distribution (mean, fraction >=10, max) in the true alignment,
and for the tree: root-to-tip CV (clock deviation) and height/mean-pairwise path ratio.
"""
import json, random, sys, statistics


def read_fasta(p):
    seqs, name = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def aln_stats(aln, pairs=2000, seed=1):
    names = list(aln); rng = random.Random(seed)
    L = len(aln[names[0]])
    gaps = sum(s.count("-") for s in aln.values())
    ung = [len(s) - s.count("-") for s in aln.values()]
    pd = []
    for _ in range(pairs):
        a, b = rng.sample(names, 2); x, y = aln[a], aln[b]
        m = d = 0
        for c1, c2 in zip(x, y):
            if c1 != "-" and c2 != "-":
                m += 1; d += c1 != c2
        if m: pd.append(d / m)
    runs = []
    for s in aln.values():
        r = 0
        for c in s:
            if c == "-": r += 1
            elif r: runs.append(r); r = 0
        if r: runs.append(r)
    return {"n": len(names), "aln_len": L, "mean_seq_len": statistics.mean(ung),
            "gap_frac": gaps / (L * len(names)), "mean_p": statistics.mean(pd), "max_p": max(pd),
            "gap_runs_per_seq": len(runs) / len(names),
            "gap_run_mean": statistics.mean(runs) if runs else 0,
            "gap_run_ge10": sum(r >= 10 for r in runs) / max(1, len(runs)),
            "gap_run_max": max(runs) if runs else 0}


def tree_stats(path):
    import dendropy
    t = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True)
    for e in t.preorder_edge_iter():
        if e.length is None: e.length = 0.0
    t.calc_node_root_distances()
    rt = [l.root_distance for l in t.leaf_node_iter()]
    tl = sum(e.length for e in t.preorder_edge_iter())
    return {"rtt_mean": statistics.mean(rt), "rtt_cv": statistics.pstdev(rt) / statistics.mean(rt),
            "rtt_max_over_min": max(rt) / min(rt), "tree_len": tl,
            "colless_like_depth": statistics.mean(len(list(l.ancestor_iter())) for l in t.leaf_node_iter())}


if __name__ == "__main__":
    out = aln_stats(read_fasta(sys.argv[1]))
    if len(sys.argv) > 2:
        out.update(tree_stats(sys.argv[2]))
    print(json.dumps(out))
