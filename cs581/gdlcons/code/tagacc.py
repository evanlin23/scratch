"""How well does ASTRAL-Pro3's own rooting + tagging recover the true root and D/S tags
on true GDL gene trees?  usage: python tagacc.py SETTING REP NFAM"""
import json, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core
from phylo import parse_newick

setting, rep, nfam = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
d = os.path.join("/opt/runs/gdlcons/data", setting, "rep%d" % rep)
sp = parse_newick(open(os.path.join(d, "species.nwk")).read())
species = sorted(sp.label[v] for v in sp.leaves())
spi = {s: i for i, s in enumerate(species)}
lines = [l for l in open(os.path.join(d, "fams.nwk")) if ";" in l][:nfam]
base = [core.true_rt(l, spi) for l in lines]
with tempfile.TemporaryDirectory() as td:
    gi, mp, out = (os.path.join(td, x) for x in ("g", "m", "o"))
    open(gi, "w").write("\n".join(rt.newick() for rt, _ in base) + "\n")
    core.write_mapping(mp, [rt for rt, _ in base])
    subprocess.run([core.APRO, "-T", "-a", mp, "-o", out, gi], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    est = [l for l in open(out) if ";" in l]
root_ok = n_int = mis = d2s = s2d = ndup = 0
ovl_mis = 0
for (rt, tg), e in zip(base, est):
    et = parse_newick(e)
    ctrue = core.clusters(rt)
    cest = core.clusters(et)
    rootset = lambda t: frozenset(frozenset(next(k for k, v in core.clusters(t).items() if v == c) if t.children[c] else [t.label[c]]) for c in t.children[t.root])
    root_ok += rootset(rt) == rootset(et)
    ov = core.overlap_tags(rt, spi)
    for cl, v in ctrue.items():
        n_int += 1
        ndup += tg[v]
        ovl_mis += ov[v] != tg[v]
        if cl in cest:
            u = cest[cl]
            et_tag = et.label[u] == "D"
            if et_tag != tg[v]:
                mis += 1
                if tg[v]: d2s += 1
                else: s2d += 1
print(json.dumps({"setting": setting, "rep": rep, "nfam": len(base), "root_correct": root_ok / len(base),
                  "internal_nodes": n_int, "true_dup_frac": ndup / n_int,
                  "own_mistag_frac_on_shared_clusters": mis / n_int, "own_D_as_S": d2s, "own_S_as_D": s2d,
                  "overlap_rule_on_true_root_mistag_frac": ovl_mis / n_int}))
