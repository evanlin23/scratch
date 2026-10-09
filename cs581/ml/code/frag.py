"""Pilot: fragment-aware two-phase ML ("backbone first") for alignments with fragmentary sequences.

    python frag.py DATASET REP [BACKBONE_METHOD]   (BACKBONE_METHOD: fasttree (default) | raxmlng)

1. backbone = full-length sequences (ungapped length >= 50% of the median);
2. tree T_b on the backbone alignment with BACKBONE_METHOD (GTR+G);
3. "frag_constr": RAxML-NG search (GTR+G, one parsimony start) on ALL sequences with T_b as a
   non-comprehensive topological constraint, so only the fragments are placed freely;
4. "frag_polish": unconstrained RAxML-NG search started from the step-3 tree.
Rows (FN vs true tree, cumulative wall/CPU seconds) are appended to results/frag.jsonl.
"""
import json
import os
import statistics
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402
import treeerr  # noqa: E402

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "frag.jsonl")


def main(ds, rep, bb_method="fasttree"):
    d = os.path.join(rt.MLDATA, ds, "R%s" % rep)
    true_tree = os.path.join(d, "true_tree.tre")
    names, seqs = rt.read_fasta(os.path.join(d, "true_align.fasta"))
    lens = [len(s.replace("-", "")) for s in seqs]
    med = statistics.median(lens)
    bb = [i for i, l in enumerate(lens) if l >= 0.5 * med]
    work = tempfile.mkdtemp(prefix="frag_%s_%s_" % (ds, rep))
    full = os.path.join(d, "true_align.clean.fasta")
    if not os.path.exists(full):
        rt.clean_alignment(os.path.join(d, "true_align.fasta"), full)
    bba = os.path.join(work, "backbone.fasta")
    rt.write_fasta(bba + ".raw", [names[i] for i in bb], [seqs[i] for i in bb])
    rt.clean_alignment(bba + ".raw", bba)
    tb = os.path.join(d, "trees", "backbone.%s.tre" % bb_method)
    os.makedirs(os.path.dirname(tb), exist_ok=True)
    w0 = os.path.join(work, "bb")
    os.makedirs(w0)
    sec0, _ = rt.estimate(bb_method, bba, tb, w0)
    cpu = rt.CPU.seconds
    base = {"dataset": ds, "rep": rep, "aln": "true_align", "backbone": bb_method, "n_backbone": len(bb),
            "backbone_fn": treeerr.error(true_tree, tb)["fn_rate"]}
    rows = []
    t1 = os.path.join(d, "trees", "true_align.frag_constr_%s.tre" % bb_method)
    w1 = os.path.join(work, "constr")
    os.makedirs(w1)
    sec1, lnl1 = rt.estimate("raxmlng", full, t1, w1, extra=["--tree-constraint", tb])
    cpu += rt.CPU.seconds
    e = treeerr.error(true_tree, t1)
    rows.append({**base, "method": "frag_constr_" + bb_method, "seconds": round(sec0 + sec1, 1),
                 "cpu_seconds": round(cpu, 1), "lnl_tool": lnl1, "fn_rate": e["fn_rate"], "rf_rate": e["rf_rate"]})
    t2 = os.path.join(d, "trees", "true_align.frag_polish_%s.tre" % bb_method)
    w2 = os.path.join(work, "polish")
    os.makedirs(w2)
    sec2, lnl2 = rt.estimate("raxmlng_ft", full, t2, w2, start_tree=t1)
    cpu += rt.CPU.seconds
    e = treeerr.error(true_tree, t2)
    rows.append({**base, "method": "frag_polish_" + bb_method, "seconds": round(sec0 + sec1 + sec2, 1),
                 "cpu_seconds": round(cpu, 1), "lnl_tool": lnl2, "fn_rate": e["fn_rate"], "rf_rate": e["rf_rate"]})
    with open(RESULTS, "a") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
            print(json.dumps(r), flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
