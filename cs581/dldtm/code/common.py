"""Shared helpers: dataset paths, FASTA I/O, subset-tree estimators, centroid decomposition."""
import os
import re
import subprocess
import tempfile

from phylo import read_tree, parse_newick

DATA = os.environ.get("DATA", "/opt/data/Datasets")
IQ = os.environ.get("IQTREE", "/opt/mm/root/envs/bio/bin/iqtree3")
FASTTREE = os.environ.get("FASTTREE", "/usr/bin/FastTree")
GTM = os.environ.get("GTM_SRC", "/opt/gtm/GTM_src")
PF_REPO = os.environ.get("PF_REPO", "/opt/src/Phyloformer")
FASTME = f"{PF_REPO}/bin/bin_linux/fastme"


def paths(cond, rep):
    """cond like 1000M2 or RNASim1000; rep like R0. Returns (true_aln, true_tree)."""
    if cond.startswith("RNASim"):
        d = f"{DATA}/RNASim/{cond[6:]}/{rep}"
        return f"{d}/true_align.txt", f"{d}/true_tree.tre"
    d = f"{DATA}/ROSE/{cond}/{rep}"
    return f"{d}/rose.aln.true.fasta", f"{d}/rose.tt"


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name and line:
            seqs[name].append(line)
    return {k: "".join(v).upper().replace("U", "T") for k, v in seqs.items()}


def write_fasta(seqs, names, path, strip_gap_cols=True):
    rows = [seqs[n] for n in names]
    if strip_gap_cols:
        keep = [i for i in range(len(rows[0])) if any(r[i] != "-" for r in rows)]
        rows = ["".join(r[i] for i in keep) for r in rows]
    with open(path, "w") as f:
        for n, r in zip(names, rows):
            f.write(f">{n}\n{r}\n")


def write_phylip_aln(fasta, path):
    s = read_fasta(fasta)
    with open(path, "w") as f:
        f.write(f"{len(s)} {len(next(iter(s.values())))}\n")
        for i, (k, v) in enumerate(s.items()):
            f.write(f"s{i} {v}\n")
    return list(s)


def _rename_back(newick, names):
    return re.sub(r"\bs(\d+)\b(?=[:,);])", lambda m: names[int(m.group(1))], newick)


def run_fasttree(aln, out, threads=1):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(out, "w") as f:
        subprocess.run([FASTTREE, "-nt", "-gtr", "-quiet", "-nopr", aln], stdout=f, check=True,
                       stderr=subprocess.DEVNULL, env=env)


def run_iqtree(aln, out, threads=1, model="GTR+G", fast=False):
    with tempfile.TemporaryDirectory() as td:
        cmd = [IQ, "-s", aln, "-m", model, "-T", str(threads), "-pre", f"{td}/x", "-quiet", "-seed", "1"]
        if fast:
            cmd.append("--fast")
        subprocess.run(cmd, check=True, capture_output=True)
        os.replace(f"{td}/x.treefile", out)


def run_fastme_seq(aln, out, method="N", dist="J", spr=False):
    """FastME on sequences: method N=NJ, B=BME; dist J=JC69, K=K2P, T=TN93, L=LogDet."""
    with tempfile.TemporaryDirectory() as td:
        names = write_phylip_aln(aln, f"{td}/a.phy")
        cmd = [FASTME, "-i", f"{td}/a.phy", "-d" + dist, "-m", method, "-o", f"{td}/t.nwk", "-T", "1"]
        if spr:
            cmd += ["--nni", "--spr"]
        subprocess.run(cmd, check=True, capture_output=True)
        s = open(f"{td}/t.nwk").read().strip()
    with open(out, "w") as f:
        f.write(_rename_back(s, names) + "\n")


def centroid_decomp(t, maxsize):
    """Recursive centroid-edge decomposition of phylo.Tree t -> list of leaf-name lists
    (same as cs581/gtm/code/sim.py, the Park et al. 2021 decomposition)."""
    names = [t.label[v] for v in t.leaves()]
    if len(names) <= maxsize:
        return [names]
    root = t.leaves()[0]
    par = {root: None}
    order = []
    st = [root]
    while st:
        v = st.pop()
        order.append(v)
        for u in t.adj[v]:
            if u not in par:
                par[u] = v
                st.append(u)
    cnt = {}
    for v in reversed(order):
        cnt[v] = (1 if v in t.label else 0) + sum(cnt[u] for u in t.adj[v] if u != par[v])
    N = len(names)
    best = min((v for v in order if par[v] is not None), key=lambda v: abs(N - 2 * cnt[v]))
    side = set()
    st = [best]
    while st:
        v = st.pop()
        if v in t.label:
            side.add(t.label[v])
        for u in t.adj[v]:
            if u != par[v]:
                st.append(u)
    a = t.copy().restrict(side)
    b = t.copy().restrict(set(names) - side)
    return centroid_decomp(a, maxsize) + centroid_decomp(b, maxsize)
