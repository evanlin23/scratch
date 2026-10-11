"""Build rogue-taxon instances on top of a ROSE replicate (true alignment and tree stay exact).

    python simulate.py inject ROSE_REP_DIR OUT_DIR --k 50 --length 1.5 --seed 1
    python simulate.py longbranch ROSE_REP_DIR OUT_DIR --k 50

inject: k rogue leaves, each evolved with AliSim (HKY+G4 + indels) along a pendant branch of
`length` substitutions/site from a random node (leaf or internal) of the ROSE model tree,
starting from that node's true sequence (rose.aln.true.internal.fasta). Rogue residues that
AliSim keeps homologous to the anchor go into the anchor's column; rogue insertions get
their own columns. The ROSE leaves' true alignment is unchanged.

longbranch: no new sequences; the rogues are the k ROSE leaves with the longest terminal
branches in the model tree.

OUT_DIR gets: true.fasta (all taxa), unaligned.fasta, true.tre, rogues.txt, info.json.
"""

import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))
sys.path.insert(0, HERE)
from gcmx import fasta  # noqa: E402
import trees  # noqa: E402

ALISIM = "/opt/mm/root/envs/bio/bin/iqtree3"


def drop_allgap(aln):
    names = list(aln)
    L = len(aln[names[0]])
    keep = [i for i in range(L) if any(aln[n][i] != "-" for n in names)]
    return {n: "".join(aln[n][i] for i in keep) for n in names}


def evolve(anchor_seq, length, seed, tmp, indel, model):
    """Evolve one sequence from anchor_seq along a branch; return (anchor_row, rogue_row) pairwise true alignment."""
    with open(os.path.join(tmp, "root.fa"), "w") as f:
        f.write(">ANC\n" + anchor_seq + "\n")
    with open(os.path.join(tmp, "t.nwk"), "w") as f:
        f.write("(COPY:0.0000001,ROGUE:%g);\n" % length)
    out = os.path.join(tmp, "sim")
    cmd = [ALISIM, "--alisim", out, "-t", os.path.join(tmp, "t.nwk"), "-m", model,
           "--root-seq", os.path.join(tmp, "root.fa") + ",ANC", "--indel", indel,
           "--indel-size", "GEO{3},GEO{3}", "-af", "fasta", "--seed", str(seed), "--no-unaligned", "-redo"]
    subprocess.run(cmd, check=True, capture_output=True, cwd=tmp)
    a = fasta.read(out + ".fa")
    return a["COPY"].upper(), a["ROGUE"].upper()


def inject(rep, out, k, length, seed, indel, model):
    rng = random.Random(seed)
    internal = fasta.upper(fasta.read(os.path.join(rep, "rose.aln.true.internal.fasta")))
    leaves = fasta.upper(fasta.read(os.path.join(rep, "rose.aln.true.fasta")))
    tree = trees.read(os.path.join(rep, "rose.mt.internal"))
    nodes = [n for n in tree.postorder() if n.parent is not None and n.name in internal]
    anchors = [rng.choice(nodes) for _ in range(k)]
    cols = len(next(iter(internal.values())))
    # per original column: residues placed there by rogues; and insertions after it (per rogue)
    rows = {}
    with tempfile.TemporaryDirectory() as tmp:
        for i, a in enumerate(anchors):
            name = "ROGUE%03d" % i
            arow = internal[a.name]
            colpos = [c for c, ch in enumerate(arow) if ch != "-"]  # column of each anchor residue
            copy, rog = evolve(arow.replace("-", ""), length, seed * 1000 + i, tmp, indel, model)
            assert copy.replace("-", "") == arow.replace("-", ""), "anchor copy changed"
            placed = {}          # column -> residue
            ins = {}             # column (insertion after it; -1 = before all) -> list of residues
            r = -1               # index of last anchor residue seen
            for ca, cr in zip(copy, rog):
                if ca != "-":
                    r += 1
                    if cr != "-":
                        placed[colpos[r]] = cr
                elif cr != "-":
                    ins.setdefault(colpos[r] if r >= 0 else -1, []).append(cr)
            rows[name] = (placed, ins, a.name)
    # assemble: original column c, then each rogue's insertion block after c in its own columns
    names = list(leaves)
    out_rows = {n: [] for n in names + list(rows)}
    order = sorted(rows)
    def emit_ins(c):
        for rn in order:
            block = rows[rn][1].get(c)
            if not block:
                continue
            for n in out_rows:
                out_rows[n].append("".join(block) if n == rn else "-" * len(block))
    emit_ins(-1)
    for c in range(cols):
        for n in names:
            out_rows[n].append(leaves[n][c])
        for rn in order:
            out_rows[rn].append(rows[rn][0].get(c, "-"))
        emit_ins(c)
    aln = drop_allgap({n: "".join(v) for n, v in out_rows.items()})
    assert all(aln[n].replace("-", "") == leaves[n].replace("-", "") for n in names)
    assert fasta.restrict(aln, names) == drop_allgap(leaves)
    # true tree: rogue attached at its anchor (sister of a leaf anchor, or extra child of an internal one)
    for rn in order:
        anchor = next(n for n in tree.postorder() if n.name == rows[rn][2])
        if anchor.children:
            anchor.add(trees.Node(rn, length))
        else:
            p = anchor.parent
            mid = trees.Node(None, anchor.length)
            p.children[p.children.index(anchor)] = mid
            mid.parent = p
            anchor.length = 0.0
            mid.add(anchor)
            mid.add(trees.Node(rn, length))
    for n in tree.postorder():
        if n.children:
            n.name = None
    save(out, aln, tree, order, {"mode": "inject", "k": k, "length": length, "seed": seed, "indel": indel,
                                 "model": model, "anchors": {rn: rows[rn][2] for rn in order},
                                 "anchor_is_leaf": sum(rows[rn][2] in leaves for rn in order)})


def longbranch(rep, out, k):
    leaves = fasta.upper(fasta.read(os.path.join(rep, "rose.aln.true.fasta")))
    tree = trees.read(os.path.join(rep, "rose.mt.internal"))
    term = sorted(((l.length, l.name) for l in tree.leaves()), reverse=True)
    rogues = sorted(n for _, n in term[:k])
    for n in tree.postorder():
        if n.children:
            n.name = None
    save(out, leaves, tree, rogues, {"mode": "longbranch", "k": k, "min_terminal": term[k - 1][0],
                                     "median_terminal": term[len(term) // 2][0]})


def save(out, aln, tree, rogues, info):
    os.makedirs(out, exist_ok=True)
    fasta.write(aln, os.path.join(out, "true.fasta"))
    fasta.write(fasta.ungap(aln), os.path.join(out, "unaligned.fasta"))
    with open(os.path.join(out, "true.tre"), "w") as f:
        f.write(trees.write(tree) + "\n")
    with open(os.path.join(out, "rogues.txt"), "w") as f:
        f.write("\n".join(rogues) + "\n")
    info.update({"ntaxa": len(aln), "ncols": len(next(iter(aln.values())))})
    json.dump(info, open(os.path.join(out, "info.json"), "w"), indent=1)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=("inject", "longbranch"))
    p.add_argument("rep")
    p.add_argument("out")
    p.add_argument("--k", type=int, default=50)
    p.add_argument("--length", type=float, default=1.5)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--indel", default="0.05,0.05")
    p.add_argument("--model", default="K80{2}+G4{1}")
    a = p.parse_args()
    if a.mode == "inject":
        inject(a.rep, a.out, a.k, a.length, a.seed, a.indel, a.model)
    else:
        longbranch(a.rep, a.out, a.k)
