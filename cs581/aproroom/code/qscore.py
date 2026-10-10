"""Is ASTRAL-Pro3's error a search failure or an objective failure?
For each input, score under ASTRAL-Pro3's objective (-C -c): the true species tree, ASTRAL-Pro3's own tree and the
ASTRID-Pro tree.  Usage: python qscore.py OUT.jsonl  (FastMulRFS 100 genes, est25 and est100, reps 01-05)"""
import json, os, re, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import mapping, run_method, rf_shared  # noqa: E402
from phylo import parse_newick  # noqa: E402
AP = "/opt/mm/root/envs/gdl/bin/astral-pro3"
out = sys.argv[1]
done = set()
if os.path.exists(out):
    done = {(r["cond"], r["rep"], r["level"]) for r in map(json.loads, open(out))}


def score(genes, mp, tree, td):
    t = os.path.join(td, "c.tre")
    open(t, "w").write(tree.strip() + "\n")
    p = subprocess.run([AP, "-C", "-c", t, "-a", mp, "-i", genes, "-o", os.path.join(td, "s.tre")],
                       capture_output=True, text=True)
    m = re.findall(r"Score: ([0-9.e+-]+)", p.stderr)
    return float(m[-1]) if m else None


for rep in ["%02d" % i for i in range(1, 6)]:
    for dl in ("0.0000000001", "0.0000000002", "0.0000000005"):
        for ps in ("10000000", "50000000"):
            for sq in (25, 100):
                c = f"dl{dl[-1]}e-10_ps{ps[0]}e7"
                if (c, rep, f"est{sq}") in done:
                    continue
                d = f"/opt/data/fmrfs/ntaxa-100.dlrate-{dl}.psize-{ps}/{rep}"
                with tempfile.TemporaryDirectory(dir="/opt/tmp") as td:
                    g = os.path.join(td, "genes.trees")
                    with open(g, "w") as f:
                        f.writelines([l for l in open(f"{d}/g_trees-raxml-sqlen-{sq}.trees") if ";" in l][:100])
                    mp = os.path.join(td, "map.txt")
                    mapping(g, mp, "simphy")
                    true = open(f"{d}/s_tree.trees").read().strip()
                    tt = parse_newick(true.replace("[&R]", ""))
                    rec = {"cond": c, "rep": rep, "level": f"est{sq}"}
                    for m in ("astral-pro3", "astrid-pro"):
                        nwk = run_method(m, g, td, "simphy")
                        rec[m + "_FN"] = rf_shared(parse_newick(nwk), tt)[0]
                        rec[m + "_score"] = score(g, mp, nwk, td)
                    rec["true_score"] = score(g, mp, true.replace("[&R] ", ""), td)
                    open(out, "a").write(json.dumps(rec) + "\n")
                    print(json.dumps(rec), flush=True)
