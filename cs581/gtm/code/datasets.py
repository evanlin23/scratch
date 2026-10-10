"""Locate the published GTM-pipeline inputs (guide tree, subset trees, alignment)
and outputs for each condition/replicate of doi:10.13012/B2IDB-7008049_V1.
DATA root defaults to /opt/gtmdata (see fetch_data.sh)."""
import glob
import os
import re

DATA = os.environ.get("GTMDATA", "/opt/gtmdata")
R = f"{DATA}/d/RNASim1000"
C = f"{DATA}/cox/Cox1-Het"
M = f"{DATA}/m1/1000M1_HF"

REPS = {"RNASim1000": [str(i) for i in range(1, 6)],
        "Cox1-HET": ["R%02d" % i for i in range(1, 11)],
        "1000M1-HF": ["R%d" % i for i in range(5)]}


def _num(p):
    return int(re.findall(r"(\d+)\.(?:treefile|out|fas)$", p)[0])


def inputs(cond, rep, guide):
    """guide in {'FT','IQ'}. Returns dict(guide=, subsets=[...], aln=[...],
    true=, published_gtm=)."""
    if cond == "RNASim1000":
        d = f"{R}/CreateConstraintTrees/500/{'FastTree' if guide == 'FT' else 'IQTree2'}/{rep}/output"
        g = f"{d}/fasttree.out" if guide == "FT" else f"{d}/iqtree-full.treefile"
        subs = glob.glob(f"{d}/sequence_partition_*.treefile")
        aln = [s.replace(".treefile", ".out") for s in subs]
        true = f"{R}/{rep}/model/true.tt"
        pub = f"{R}/GTM/{'fasttree_fasttree' if guide == 'FT' else 'IQTree2'}/500/{rep}/branch_length."
    elif cond == "Cox1-HET":
        d = f"{C}/ConstraintTrees/{rep}/{'fasttree' if guide == 'FT' else 'iqtree+GTR+G'}-centro-500"
        g = f"{d}/fasttree.tre" if guide == "FT" else f"{d}/iqtree-result.treefile"
        subs = [p for p in glob.glob(f"{d}/centro-500-*.treefile") if re.search(r"-\d+\.treefile$", p)]
        aln = [s.replace(".treefile", ".fas") for s in subs]
        true = f"{DATA}/cox/Cox1-HET/{rep}/true-tree.tre"
        pub = f"{C}/GTM/500/{'fasttree_fasttree' if guide == 'FT' else 'iqtree'}/{rep}/branch_length."
    elif cond == "1000M1-HF":
        d = f"{M}/CreateConstraintTrees/{'FastTree' if guide == 'FT' else 'IQTree2'}/500/{rep}/output"
        g = f"{d}/fasttree.out" if guide == "FT" else f"{d}/iqtree-full.treefile"
        subs = glob.glob(f"{d}/sequence_partition_*.treefile")
        aln = [s.replace(".treefile", ".out") for s in subs]
        true = f"{DATA}/m1_true/{rep}/rose.tt"
        pub = f"{M}/GTM/500/{'fasttree_fasttree' if guide == 'FT' else 'iqtree'}/{rep}/branch_length."
    else:
        raise ValueError(cond)
    subs.sort(key=_num)
    aln.sort(key=_num)
    return dict(guide=g, subsets=subs, aln=aln, true=true, published_gtm=pub)


def read_fasta(paths):
    seqs = {}
    for p in paths:
        name = None
        with open(p) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    seqs[name] = []
                elif name:
                    seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}
