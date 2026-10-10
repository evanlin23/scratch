"""TWILIGHT iterative mode, re-implemented from its Snakemake workflow (TWILIGHT repo, workflow/).

    python3 twilight_iter.py UNALIGNED OUT THREADS WORKDIR [ITER]

TWILIGHT 0.2.3's CLI needs a guide tree. The official workflow (config.yaml defaults:
3 iterations, rgc 0.95, gap-open -50, gap-extend -5) builds tree_iter0 with DIPPER,
MAFFT --treeout or MAFFT PartTree, aligns with TWILIGHT, then re-estimates the tree with
FastTree (-gtr -nt -fastest on the alignment with >=95%-gap columns masked) and realigns.
DIPPER is not installed here, so the initial tree is MAFFT PartTree (the workflow's
"parttree" option) and later trees use FastTree (its "fasttree" option).
Polytomies are resolved randomly with DendroPy (the workflow uses ete3).
"""

import os
import subprocess
import sys

BIO = "/opt/mm/root/envs/bio/bin"
WF = "/opt/tools/twilight_wf"  # mafft2nwk.py from the TWILIGHT repo


def sh(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def read(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line:
            seqs[name].append(line)
    return {n: "".join(s) for n, s in seqs.items()}


def mask(msa, out, thr=0.95):
    seqs = read(msa)
    rows = list(seqs.values())
    n = len(rows)
    keep = [i for i in range(len(rows[0])) if sum(r[i] == "-" for r in rows) / n < thr]
    with open(out, "w") as f:
        for name, s in seqs.items():
            f.write(">{}\n{}\n".format(name, "".join(s[i] for i in keep)))


def resolve(intree, outtree):
    code = ("import dendropy,sys;t=dendropy.Tree.get(path=sys.argv[1],schema='newick',preserve_underscores=True);"
            "t.resolve_polytomies();t.write(path=sys.argv[2],schema='newick',suppress_rooting=True)")
    sh([os.path.join(BIO, "python"), "-c", code, intree, outtree])


def main():
    unal, out, threads, work = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    iters = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    os.makedirs(work, exist_ok=True)
    allseq = "".join(read(unal).values()).upper()
    nt = sum(allseq.count(c) for c in "ACGTUN") / max(1, len(allseq))
    model = ["-gtr", "-nt"] if nt > 0.9 else []  # 16S has IUPAC ambiguity codes: use the ACGTUN fraction
    seq = os.path.join(work, "seqs.fa")
    with open(seq, "w") as f:  # mafft2nwk maps PartTree's 1-based indices back to names in input order
        for n, s in read(unal).items():
            f.write(">{}\n{}\n".format(n, s))
    with open(os.path.join(work, "mafft.out"), "w") as mo:
        sh(["mafft", "--retree", "0", "--treeout", "--parttree", "--reorder", "--quiet", "--thread", threads, seq],
           stdout=mo, cwd=work)
    tree = os.path.join(work, "tree_iter0.nwk")
    sh(["python3", os.path.join(WF, "mafft2nwk.py"), seq + ".tree", seq, tree, "--parttree"])
    for it in range(1, iters + 1):
        msa = os.path.join(work, "msa_iter{}.fa".format(it))
        sh([os.path.join(BIO, "twilight"), "-i", seq, "-t", tree, "-o", msa, "-C", threads, "-r", "0.95",
            "--gap-open", "-50", "--gap-extend", "-5", "--overwrite"], cwd=work)
        if it < iters:
            masked = os.path.join(work, "mask.fa")
            mask(msa, masked)
            raw = os.path.join(work, "ft.nwk")
            with open(raw, "w") as fo:
                sh([os.path.join(BIO, "fasttreeMP"), *model, "-fastest", "-quiet", masked], stdout=fo,
                   env={**os.environ, "OMP_NUM_THREADS": threads})
            tree = os.path.join(work, "tree_iter{}.nwk".format(it))
            resolve(raw, tree)
    os.replace(msa, out)


if __name__ == "__main__":
    main()
