"""Validate apro.py against the ASTRAL-Pro3 binary on family-R gene trees: total quartet
support of the three topologies (stock binary, own rooting) vs our re-implementation."""
import os, random, subprocess, sys, tempfile
from collections import defaultdict
from gdlsim import to_newick
from famR_sim import family
from apro import family_scores

APRO = "/opt/mm/root/envs/gdl/bin/astral-pro3"
SPT = "(((A,B),C),D);"


def split_of(label):
    a, b = label.split("#")
    side = lambda x: "".join(sorted(y.strip("{}") for y in x.split("|")))
    k = "|".join(sorted([side(a), side(b)]))
    return k


def binary_scores(newicks, workdir=None, binary=APRO, env_extra=None):
    with tempfile.TemporaryDirectory(dir=workdir) as td:
        gi, mp, st, out = (os.path.join(td, x) for x in ("g.nwk", "map.txt", "sp.nwk", "s.nwk"))
        open(gi, "w").write("\n".join(newicks) + "\n")
        open(st, "w").write(SPT + "\n")
        labs = set()
        for nw in newicks:
            for tok in nw.replace("(", " ").replace(")", " ").replace(",", " ").replace(";", " ").split():
                tok = tok.rstrip("DS") if False else tok
                if "_" in tok:
                    labs.add(tok.split(":")[0])
        open(mp, "w").write("".join("%s %s\n" % (l, l.rsplit("_", 1)[0]) for l in sorted(labs)))
        env = dict(os.environ); env.pop("APRO_FIXED", None); env.update(env_extra or {})
        subprocess.run([binary, "-a", mp, "-c", st, "-C", "-u", "3", "-o", out, gi], check=True, cwd=td,
                       env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(td, "freqQuad.csv"))]
        res = {}
        for r in rows:
            k = split_of(r[2])
            if "D" in k.split("|")[0] and "A" not in k.split("|")[0]:
                pass
            res[k] = float(r[4])
        return res


if __name__ == "__main__":
    a, b, c, lamT, n, seed = map(float, sys.argv[1:7])
    rng = random.Random(int(seed))
    fams = []
    while len(fams) < n:
        g = family(lamT, 1.0, a, b, c, rng, [0])
        if g is None or g.kind == "L":
            continue
        fams.append(g)
    mine = defaultdict(float)
    for g in fams:
        for k, v in family_scores(g, "own").items():
            if "|" in k:
                mine[k] += v
    nw = [to_newick(g).replace(")D", ")").replace(")S", ")") for g in fams]
    print("binary:", binary_scores(nw))
    print("apro.py own (tie-average):", dict(mine))
