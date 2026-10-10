"""Generate nested gene-family datasets (the first K families of a file are the K-family
dataset) for the error-vs-families curves.

GDL settings use the pure-GDL simulator (true root + true D/S tags, no ILS).
DLCOAL settings use SimPhy 1.0.2 (locus tree under GDL, gene tree under the
multilocus coalescent); there is no unambiguous true tagging, so they are only used
for the method atlas.

Leaves are relabelled '<species>_<k>'; species names never contain '_'.

usage: python gen_data.py SETTING REP [NFAM]
output: /opt/runs/gdlcons/data/SETTING/repREP/{species.nwk, fams.nwk}
"""
import glob
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gdlsim import simulate, yule_tree  # noqa: E402
from phylo import parse_newick  # noqa: E402

ROOT = os.environ.get("GDLCONS_DATA", "/opt/runs/gdlcons/data")
SIMPHY = "/opt/tools/SimPhy_1.0.2/bin/simphy_lnx64"
NTAX = 16

# GDL: (internal (lam, mu), terminal (lam, mu), root branch (len, lam, mu))
GDL = {
    "gdl-crit": ((1.0, 1.0), (1.0, 1.0), (0.0, 0, 0)),
    "gdl-high": ((3.0, 3.0), (3.0, 3.0), (0.0, 0, 0)),
    "gdl-adv": ((3.0, 0.5), (0.0, 3.0), (0.0, 0, 0)),       # supercritical internal, lossy tips
    "gdl-root": ((0.0, 1.0), (0.0, 2.0), (0.3, 6.0, 0.0)),  # burst of copies above the root, then loss
    "gdl50-lossy": ((1.0, 1.5), (0.0, 4.0), (0.0, 0, 0)),   # 50 taxa: families cover few species
    "gdl50-adv": ((2.0, 0.5), (0.0, 3.0), (0.0, 0, 0)),
}
# DLCOAL (SimPhy): species tree height (generations), Ne, dup rate, loss rate (per generation)
DLCOAL = {
    "dlc-mod": (2e6, 2e5, 5e-7, 5e-7),
    "dlc-ils": (2e6, 1e6, 5e-7, 5e-7),
    "dlc-high": (2e6, 2e5, 1.5e-6, 1.5e-6),
    "dlc-asym": (2e6, 2e5, 1.5e-6, 3e-7),
    "dlc50-ils": (2e6, 1e6, 1e-6, 1e-6),
}
NTAXA = {"gdl50-lossy": 50, "gdl50-adv": 50, "dlc50-ils": 50}


def gen_gdl(setting, rep, nfam, out):
    inner, term, (rlen, rl, rm) = GDL[setting]
    sp = yule_tree(NTAXA.get(setting, NTAX), seed=7000 + rep, height=1.0)
    st = parse_newick(sp)
    lam, mu = [0.0] * len(st.parent), [0.0] * len(st.parent)
    for v in range(len(st.parent)):
        if st.parent[v] < 0:
            lam[v], mu[v] = rl, rm
        elif not st.children[v]:
            lam[v], mu[v] = term
        else:
            lam[v], mu[v] = inner
    fams, tries, over = simulate(st, lam, mu, nfam, 100000 * rep + 17, min_species=4,
                                 root_len=rlen, cap=4000, max_over=nfam)
    with open(os.path.join(out, "species.nwk"), "w") as f:
        f.write(sp + "\n")
    with open(os.path.join(out, "fams.nwk.tmp"), "w") as f:
        f.write("\n".join(fams) + "\n")
    os.rename(os.path.join(out, "fams.nwk.tmp"), os.path.join(out, "fams.nwk"))
    return {"tries": tries, "overflow": over}


_LEAF = re.compile(r"([(,])(\d+)_(\d+)_(\d+)(?=[:,)])")


def gen_dlcoal(setting, rep, nfam, out):
    height, ne, lb, ld = DLCOAL[setting]
    tmp = os.path.join(out, "simphy")
    nsim = int(nfam * 1.6)  # SimPhy does not filter on #species; oversample
    if not os.path.exists(os.path.join(tmp, "1", "s_tree.trees")):
        subprocess.run([SIMPHY, "-sl", "f:%d" % NTAXA.get(setting, NTAX), "-sb", "f:0.000001", "-st", "f:%g" % height,
                        "-sp", "f:%g" % ne, "-su", "f:0.00000001", "-lb", "f:%g" % lb, "-ld", "f:%g" % ld,
                        "-rs", "1", "-rl", "f:%d" % nsim, "-rg", "1", "-o", tmp, "-cs", str(9000 + rep),
                        "-v", "0", "-oc", "1", "-om", "0", "-od", "0"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    sp = open(os.path.join(tmp, "1", "s_tree.trees")).read().strip()
    sp = re.sub(r"([(,])(\d+)(?=:)", lambda m: m.group(1) + "T" + m.group(2), sp)
    fams = []
    files = sorted(glob.glob(os.path.join(tmp, "1", "g_trees*.trees")),
                   key=lambda p: int(re.findall(r"(\d+)\.trees$", p)[0]))
    for p in files:
        g = open(p).read().strip()
        g = _LEAF.sub(lambda m: "%sT%s_%s.%s" % (m.group(1), m.group(2), m.group(3), m.group(4)), g)
        t = parse_newick(g)
        spp = set(t.label[v].split("_")[0] for v in t.leaves())
        if len(spp) >= 4:
            fams.append(g)
        if len(fams) >= nfam:
            break
    with open(os.path.join(out, "species.nwk"), "w") as f:
        f.write(sp + "\n")
    with open(os.path.join(out, "fams.nwk.tmp"), "w") as f:
        f.write("\n".join(fams) + "\n")
    os.rename(os.path.join(out, "fams.nwk.tmp"), os.path.join(out, "fams.nwk"))
    subprocess.run(["rm", "-rf", tmp])
    return {"nfam": len(fams), "simulated": nsim}


if __name__ == "__main__":
    setting, rep = sys.argv[1], int(sys.argv[2])
    nfam = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
    out = os.path.join(ROOT, setting, "rep%d" % rep)
    os.makedirs(out, exist_ok=True)
    if os.path.exists(os.path.join(out, "fams.nwk")):
        sys.exit(0)
    info = gen_gdl(setting, rep, nfam, out) if setting in GDL else gen_dlcoal(setting, rep, nfam, out)
    print(setting, rep, info, flush=True)
