#!/usr/bin/env python3
"""Run all baselines and candidate methods on one HGT+ILS replicate (Davidson et al. 2015 data).

usage: run_rep.py <model_dir_name> <rep> <ngenes> <genetype: est|true> <outroot>
       run_rep.py --files <gene_trees> <true_tree> <ngenes> <outdir>   (any other dataset)
Writes <outroot>/<model>/<rep>/<ngenes>-<genetype>/result.json (skips work already done).
"""
import json, os, re, subprocess, sys, time

DATA = "/opt/data/hgt/HGT-Davidson-2015"
BIN = "/opt/mm/root/envs/phy/bin"
QT = "/opt/runs/qtool"
WQFM = "/opt/tools/wQFM-TREE/wQFM-TREE"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from treeutil import bipartitions, rf, read_newick_list, astrid_matrix  # noqa: E402


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, **kw)


def main():
    if sys.argv[1] == "--files":
        src, true_file, ng, d = sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
        gtype = "raxml"  # no support values -> wASTRAL skipped
    else:
        model, rep, ng, gtype, outroot = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
        d = os.path.join(outroot, model, rep, f"{ng}-{gtype}")
        src = os.path.join(DATA, model, rep, "estimatedgenetre" if gtype == "est" else "truegenetrees")
        true_file = os.path.join(DATA, model, rep, "s_tree.trees")
    os.makedirs(d, exist_ok=True)
    resf = os.path.join(d, "result.json")
    res = json.load(open(resf)) if os.path.exists(resf) else {}
    # drop a root branch length / label (e.g. RAxML "...):0.0;"), which wQFM-TREE cannot parse
    genes = [re.sub(r"\)[^()]*;$", ");", t) for t in read_newick_list(src)[:ng]]
    g = os.path.join(d, "genes.tre")
    with open(g, "w") as f:
        f.write("\n".join(genes) + "\n")
    true_tree = read_newick_list(true_file)[0]
    tb = bipartitions(true_tree)

    def save():
        json.dump(res, open(resf, "w"), indent=1)

    def record(name, tree_file, sec):
        t = open(tree_file).read().strip().split("\n")[0]
        sc = sh(f"{QT} score -g {g} -t {tree_file}").stdout.split()
        fn, fp, n = rf(tb, bipartitions(t))
        res[name] = dict(score=float(sc[0]), nscore=float(sc[1]), fn=fn, fp=fp, nint=n, sec=sec)
        save()

    def run(name, cmd, out):
        if name in res:
            return
        t0 = time.time()
        sh(cmd)
        record(name, out, time.time() - t0)

    p = lambda x: os.path.join(d, x)
    # quartet score of the true species tree (diagnostic: is the MQSST objective aligned with accuracy?)
    if "truetree" not in res:
        with open(p("true.tre"), "w") as f:
            f.write(true_tree + "\n")
        record("truetree", p("true.tre"), 0.0)
    # ---------------- baselines
    run("astral4", f"{BIN}/astral4 -t 1 -u 0 -o {p('astral4.tre')} {g}", p("astral4.tre"))
    if gtype == "est":  # wASTRAL needs gene-tree support values (FastTree SH-like)
        run("wastral", f"{BIN}/wastral -t 1 -u 0 -o {p('wastral.tre')} {g}", p("wastral.tre"))
    run("treeqmc", f"{BIN}/tree-qmc -i {g} -o {p('treeqmc.tre')}", p("treeqmc.tre"))
    run("wqfm", f"cd {WQFM} && bash run.sh {g} {p('wqfm.tre')}", p("wqfm.tre"))
    if "astrid" not in res:
        t0 = time.time()
        names, M = astrid_matrix(genes)
        with open(p("astrid.phy"), "w") as f:
            f.write(f"{len(names)}\n")
            for i, nm in enumerate(names):
                f.write(nm + " " + " ".join(f"{x:.6f}" for x in M[i]) + "\n")
        sh(f"{BIN}/fastme -i {p('astrid.phy')} -o {p('astrid.tre')} -m B -s -n >/dev/null")
        record("astrid", p("astrid.tre"), time.time() - t0)
    # ---------------- (a) SPR hill-climbing on the explicit quartet table, from each start tree
    for start in ["astral4", "treeqmc", "wqfm", "astrid"]:
        run(f"ls[{start}]", f"{QT} search -g {g} -t {p(start + '.tre')} > {p('ls_' + start + '.tre')}",
            p("ls_" + start + ".tre"))
    # ---------------- (b) ASTRAL-IV DP over clusters harvested from the other heuristics
    if "harvest" not in res:
        with open(p("guide.tre"), "w") as f:
            for s in ["treeqmc", "wqfm", "astrid", "wastral"]:
                if not os.path.exists(p(s + ".tre")):
                    continue
                f.write(open(p(s + ".tre")).read().strip().split("\n")[0] + "\n")
    run("harvest", f"{BIN}/astral4 -t 1 -u 0 -g {p('guide.tre')} -o {p('harvest.tre')} {g}", p("harvest.tre"))
    # ASTRAL-III: exact DP over its default X, and over X enlarged with the harvested trees (-e)
    A3 = "java -Xmx3g -jar /opt/mm/root/envs/phy/share/astral-tree-5.7.8-1/astral.5.7.8.jar"
    run("astral3", f"{A3} -t 0 -i {g} -o {p('astral3.tre')} 2>/dev/null", p("astral3.tre"))
    if "astral3+e" not in res:  # ASTRAL-III rejects extra trees with polytomies: resolve arbitrarily
        import treeswift
        with open(p("guide3.tre"), "w") as f:
            for t in read_newick_list(p("guide.tre")):
                ts = treeswift.read_tree_newick(t)
                ts.resolve_polytomies()
                f.write(ts.newick() + "\n")
    run("astral3+e", f"{A3} -t 0 -i {g} -e {p('guide3.tre')} -o {p('astral3e.tre')} 2>/dev/null", p("astral3e.tre"))
    # ---------------- (c) HGT-aware variants
    for mode in ["capminor", "vote"]:
        run(f"{mode}", f"{QT} search -m {mode} -g {g} -t {p('astral4.tre')} > {p(mode + '.tre')}", p(mode + ".tre"))
    if "reweight" not in res:
        t0 = time.time()
        sup = [float(x) for x in sh(f"{QT} genesupport -g {g} -t {p('astral4.tre')}").stdout.split()]
        # gene weight: integer copies 0..4 proportional to (agreement / max agreement)^4
        mx = max(sup)
        cop = [round(4 * (s / mx) ** 4) for s in sup]
        with open(p("rw_genes.tre"), "w") as f:
            for gt, c in zip(genes, cop):
                f.write((gt + "\n") * c)
        sh(f"{BIN}/astral4 -t 1 -u 0 -o {p('reweight.tre')} {p('rw_genes.tre')}")
        record("reweight", p("reweight.tre"), time.time() - t0)
        res["reweight"]["copies"] = cop
        save()


if __name__ == "__main__":
    main()
