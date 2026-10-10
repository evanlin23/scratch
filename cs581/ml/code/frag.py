"""Fragment-aware two-phase ML ("backbone first") for alignments with fragmentary sequences.

    python frag.py DATASET REP [--backbone B] [--polish] [--out results/frag.jsonl] [--tau 0.5]

1. backbone = sequences with ungapped length >= TAU * median (the full-length ones);
2. backbone tree T_b with method B:
     fasttree | iqtree_fast | iqtree | raxmlng (1 parsimony start) |
     raxmlng_ft (RAxML-NG started from the FastTree backbone tree) |
     true (oracle: the true tree restricted to the backbone -- an upper bound, not a method);
3. frag_constr_<B>: RAxML-NG (GTR+G, one parsimony start) on ALL sequences with T_b as a
   non-comprehensive topological constraint, so only the fragments are placed freely;
4. --polish: frag_polish_<B>: unconstrained RAxML-NG search started from the step-3 tree.
Rows (FN vs true tree, cumulative wall/CPU seconds, backbone FN) are appended to --out.
"""
import argparse
import json
import os
import statistics
import sys
import tempfile

import dendropy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402
import treeerr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def backbone_tree(method, bba, tb, work, true_tree, labels):
    """Estimate the backbone tree; returns (wall seconds, cpu seconds)."""
    if method == "true":
        t = dendropy.Tree.get(path=true_tree, schema="newick", preserve_underscores=True)
        t.retain_taxa_with_labels(set(labels))
        t.write(path=tb, schema="newick", suppress_rooting=True)
        return 0.0, 0.0
    if method == "raxmlng_ft":
        ft = tb + ".fasttree"
        w = os.path.join(work, "bbft")
        os.makedirs(w)
        s0, _ = rt.estimate("fasttree", bba, ft, w)
        c0 = rt.CPU.seconds
        s1, _ = rt.estimate("raxmlng_ft", bba, tb, w, start_tree=ft)
        return s0 + s1, c0 + rt.CPU.seconds
    s, _ = rt.estimate(method, bba, tb, work)
    return s, rt.CPU.seconds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("rep")
    ap.add_argument("--backbone", default="fasttree")
    ap.add_argument("--polish", action="store_true")
    ap.add_argument("--tau", type=float, default=0.5)
    ap.add_argument("--aln", default="true_align", help="alignment name in the replicate dir (e.g. upp)")
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "frag.jsonl"))
    a = ap.parse_args()
    ds, rep, bbm = a.dataset, a.rep, a.backbone
    tag = bbm if a.tau == 0.5 else "%s_tau%g" % (bbm, a.tau)
    d = os.path.join(rt.MLDATA, ds, "R%s" % rep)
    true_tree = os.path.join(d, "true_tree.tre")
    names, seqs = rt.read_fasta(os.path.join(d, a.aln + ".fasta"))
    lens = [len(s.replace("-", "")) for s in seqs]
    med = statistics.median(lens)
    bb = [i for i, l in enumerate(lens) if l >= a.tau * med]
    work = tempfile.mkdtemp(prefix="frag_%s_%s_%s_" % (ds, rep, tag))
    full = os.path.join(d, a.aln + ".clean.fasta")
    if not os.path.exists(full):
        tmp = full + ".%d.tmp" % os.getpid()
        rt.clean_alignment(os.path.join(d, a.aln + ".fasta"), tmp)
        os.replace(tmp, full)
    bba = os.path.join(work, "backbone.fasta")
    rt.write_fasta(bba + ".raw", [names[i] for i in bb], [seqs[i] for i in bb])
    rt.clean_alignment(bba + ".raw", bba)
    os.makedirs(os.path.join(d, "trees"), exist_ok=True)
    tb = os.path.join(d, "trees", "%sbackbone.%s.tre" % ("" if a.aln == "true_align" else a.aln + ".", tag))
    w0 = os.path.join(work, "bb")
    os.makedirs(w0)
    sec, cpu = backbone_tree(bbm, bba, tb, w0, true_tree, [names[i] for i in bb])
    base = {"dataset": ds, "rep": rep, "aln": a.aln, "backbone": bbm, "tau": a.tau, "n_backbone": len(bb),
            "backbone_fn": treeerr.error(true_tree, tb)["fn_rate"], "backbone_cpu": round(cpu, 1)}
    rows = []
    t1 = os.path.join(d, "trees", "%s.frag_constr_%s.tre" % (a.aln, tag))
    w1 = os.path.join(work, "constr")
    os.makedirs(w1)
    s1, lnl1 = rt.estimate("raxmlng", full, t1, w1, extra=["--tree-constraint", tb])
    sec, cpu = sec + s1, cpu + rt.CPU.seconds
    e = treeerr.error(true_tree, t1)
    rows.append({**base, "method": "frag_constr_" + tag, "seconds": round(sec, 1), "cpu_seconds": round(cpu, 1),
                 "lnl_tool": lnl1, "fn_rate": e["fn_rate"], "rf_rate": e["rf_rate"], "tree": t1})
    if a.polish:
        t2 = os.path.join(d, "trees", "%s.frag_polish_%s.tre" % (a.aln, tag))
        w2 = os.path.join(work, "polish")
        os.makedirs(w2)
        s2, lnl2 = rt.estimate("raxmlng_ft", full, t2, w2, start_tree=t1)
        sec, cpu = sec + s2, cpu + rt.CPU.seconds
        e = treeerr.error(true_tree, t2)
        rows.append({**base, "method": "frag_polish_" + tag, "seconds": round(sec, 1), "cpu_seconds": round(cpu, 1),
                     "lnl_tool": lnl2, "fn_rate": e["fn_rate"], "rf_rate": e["rf_rate"], "tree": t2})
    with open(a.out, "a") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
            print(json.dumps(r), flush=True)


if __name__ == "__main__":
    main()
