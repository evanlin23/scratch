"""Pools of informative (all 4 species present) 4-taxon gene families with true D/S labels.
usage: python pools.py NAME N   -> /opt/runs/ap2/NAME.nwk  (restartable: skips if present)"""
import os, random, sys
from gdlsim import to_newick, simulate
from phylo import parse_newick
from famR_sim import family
from apro import to_unrooted

CONFIGS = {
    "R1": ("famR", (0.01, 0.01, 0.2, 3.0)),
    "R2": ("famR", (0.05, 0.05, 0.2, 3.0)),
    "cand4": ("gdl", "A=0,8,0.5;B=0,8,0.5;C=8,8,0.2;x=1,0.5,0.5;y=4,1,1;D=0,0,1"),
    "ctrl": ("gdl", "A=1,1,0.3;B=1,1,0.3;C=1,1,0.6;x=1,1,0.3;y=1,1,0.4;D=1,1,1"),
}
OUT = "/opt/runs/ap2"


def make(name, n, seed=4242):
    kind, par = CONFIGS[name]
    path = os.path.join(OUT, name + ".nwk")
    if os.path.exists(path):
        return path
    os.makedirs(OUT, exist_ok=True)
    fams = []
    if kind == "famR":
        a, b, c, lamT = par
        rng = random.Random(seed)
        while len(fams) < n:
            g = family(lamT, 1.0, a, b, c, rng, [0])
            if g is None or g.kind == "L":
                continue
            _, sp, _ = to_unrooted(g)
            if len(set(s for s in sp if s)) == 4:
                fams.append(to_newick(g))
    else:
        br = {k: tuple(float(x) for x in v.split(",")) for k, v in (kv.split("=") for kv in par.split(";"))}
        nwk = "(((A:%g,B:%g)x:%g,C:%g)y:%g,D:%g);" % (br["A"][2], br["B"][2], br["x"][2], br["C"][2], br["y"][2], br["D"][2])
        st = parse_newick(nwk)
        lam, mu = [0.0] * len(st.parent), [0.0] * len(st.parent)
        for v in range(len(st.parent)):
            if st.label[v] in br:
                lam[v], mu[v] = br[st.label[v]][0], br[st.label[v]][1]
        fams, _, _ = simulate(st, lam, mu, n, seed, min_species=4, cap=6000, max_over=n, max_tries=3000 * n)
    open(path + ".tmp", "w").write("\n".join(fams) + "\n")
    os.rename(path + ".tmp", path)
    return path


if __name__ == "__main__":
    print(make(sys.argv[1], int(sys.argv[2])))
