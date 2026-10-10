"""Species-tree error vs number of gene families, for many methods and (GDL only)
ASTRAL-Pro under rooting/tagging error models.

usage: python run_curve.py SETTING REP OUT.jsonl [GROUP]
GROUP: 'errors' (ASTRAL-Pro error models; GDL settings only), 'atlas' (methods), 'all'
Restartable: (setting, rep, nfam, method) already in OUT are skipped.
"""
import json
import os
import random
import subprocess
import sys
import tempfile
import time
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core  # noqa: E402
import methods as M  # noqa: E402
from phylo import parse_newick  # noqa: E402

DATA = os.environ.get("GDLCONS_DATA", "/opt/runs/gdlcons/data")
KS = [int(x) for x in os.environ.get("KS", "50,200,1000,5000,20000").split(",")]
PY = "/opt/mm/root/envs/gdl/bin/python"
PY27 = "/opt/mm/root/envs/py27/bin/python"
DISCO = "/opt/tools/DISCO/disco.py"
FMPRE = "/opt/tools/fastmulrfs/python-tools/preprocess_multrees_v3.py"
FASTRFS = "/opt/tools/fastmulrfs/external/FastRFS/build/FastRFS"
STAG = "/opt/tools/STAG/stag/stag.py"
TIMEOUT = int(os.environ.get("METHOD_TIMEOUT", "2400"))

ERROR_MODELS = [("true", 0), ("ovl", 0), ("rovl", 0.1), ("rovl", 0.3), ("rovl", 1.0),
                ("flip", 0.05), ("flip", 0.15), ("flip", 0.3), ("d2s", 0.5), ("d2s", 1.0),
                ("s2d", 0.5)]
ATLAS = ["astral-pro", "astrid-multi", "astrid-disco", "astral-disco", "fastmulrfs", "stag"]


def species_of_tree(t):
    return sorted(set(core.sp_of(t.label[v]) for v in t.leaves()))


def run_disco(plain_file, out_file):
    subprocess.run([PY, DISCO, "-i", plain_file, "-o", out_file, "-d", "_", "-m", "4"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=TIMEOUT)
    return [l.strip() for l in open(out_file) if ";" in l]


def astrid_single(trees_nwk, species):
    """ASTRID (average internode distance + FastME) on single-copy trees."""
    spi = {s: i for i, s in enumerate(species)}
    per = []
    for nw in trees_nwk:
        t = parse_newick(nw)
        tg = [False] * len(t.parent)
        ls = [spi[t.label[v]] if not t.children[v] else -1 for v in range(len(t.parent))]
        per.append(M.gene_distances(t, tg, ls, len(species), mode="multi"))
    D, _ = M.average_matrix(per)
    return M.fastme_tree(D, species)


def method_tree(method, fams_plain, rts, species, td):
    gi = os.path.join(td, "g.nwk")
    with open(gi, "w") as f:
        f.write("\n".join(fams_plain) + "\n")
    if method == "astral-pro":
        return core.run_apro(fams_plain, rts, fixed=False)
    if method == "astrid-multi":
        spi = {s: i for i, s in enumerate(species)}
        per = []
        for rt in rts:
            ls = [spi[core.sp_of(rt.label[v])] if not rt.children[v] else -1 for v in range(len(rt.parent))]
            per.append(M.gene_distances(rt, [False] * len(rt.parent), ls, len(species), mode="multi"))
        D, _ = M.average_matrix(per)
        return M.fastme_tree(D, species)
    if method in ("astrid-disco", "astral-disco"):
        dfile = os.path.join(td, "disco.nwk")
        if not os.path.exists(dfile):
            run_disco(gi, dfile)
        sc = [l.strip() for l in open(dfile) if ";" in l]
        if method == "astrid-disco":
            return astrid_single(sc, species)
        out = os.path.join(td, "ad.nwk")
        subprocess.run([core.ASTRAL, "-i", dfile, "-o", out, "-t", "1"], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=TIMEOUT)
        return parse_newick(open(out).read())
    if method == "fastmulrfs":
        # FastMulRFS wants leaves labelled by species
        sp_lab = os.path.join(td, "g_sp.nwk")
        with open(sp_lab, "w") as f:
            for rt in rts:
                t2 = parse_newick(rt.newick())
                for v in t2.leaves():
                    t2.label[v] = core.sp_of(t2.label[v])
                f.write(t2.newick() + "\n")
        pre = os.path.join(td, "g_fr.nwk")
        subprocess.run([PY, FMPRE, "-i", sp_lab, "-o", pre], check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=TIMEOUT)
        out = os.path.join(td, "fm.tree")
        subprocess.run([FASTRFS, "-i", pre, "-o", out], check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=TIMEOUT, cwd=td)
        return parse_newick(open(out + ".single").read())
    if method == "stag":
        gdir = os.path.join(td, "stag_trees")
        os.makedirs(gdir, exist_ok=True)
        n_all = 0
        for i, nw in enumerate(fams_plain):
            t = parse_newick(nw)
            if len(species_of_tree(t)) == len(species):
                n_all += 1
                with open(os.path.join(gdir, "OG%07d_tree.txt" % i), "w") as f:
                    f.write(nw + "\n")
        if n_all == 0:
            raise RuntimeError("no family contains all species")
        smap = os.path.join(td, "smap.txt")
        with open(smap, "w") as f:
            for s in species:
                f.write("%s_* %s\n" % (s, s))
        env = dict(os.environ, PATH="/opt/mm/root/envs/gdl/bin:" + os.environ["PATH"])
        subprocess.run([PY27, STAG, smap, gdir], check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=TIMEOUT, cwd=td, env=env)
        res = [os.path.join(dp, f) for dp, _, fs in os.walk(td) for f in fs if f == "SpeciesTree.tre"]
        return parse_newick(open(res[0]).read())
    raise ValueError(method)


def main():
    setting, rep, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    group = sys.argv[4] if len(sys.argv) > 4 else "all"
    d = os.path.join(DATA, setting, "rep%d" % rep)
    sp_nwk = open(os.path.join(d, "species.nwk")).read().strip()
    lines = [l.strip() for l in open(os.path.join(d, "fams.nwk")) if ";" in l]
    true_sp = parse_newick(sp_nwk)
    species = sorted(true_sp.label[v] for v in true_sp.leaves())
    spi = {s: i for i, s in enumerate(species)}
    done = set()
    if os.path.exists(out):
        for l in open(out):
            r = json.loads(l)
            done.add((r["setting"], r["rep"], r["nfam"], r["method"]))
    gdl = setting.startswith("gdl")
    for k in KS:
        if k > len(lines):
            continue
        todo = []
        if group in ("errors", "all") and gdl:
            todo += ["apro-fixed:%s:%g" % m for m in ERROR_MODELS]
        if group in ("atlas", "all"):
            todo += ATLAS
        todo = [m for m in todo if (setting, rep, k, m) not in done]
        if not todo:
            continue
        sub = lines[:k]
        if gdl:
            base = [core.true_rt(nw, spi) for nw in sub]
        else:
            base = [(parse_newick(nw), None) for nw in sub]
        rts = [rt for rt, _ in base]
        plain = [rt.newick() for rt in rts]
        with tempfile.TemporaryDirectory(dir="/opt/runs/gdlcons/tmp") as td:
            for m in todo:
                t0 = time.time()
                rec = {"setting": setting, "rep": rep, "nfam": k, "method": m}
                try:
                    if m.startswith("apro-fixed:"):
                        _, model, par = m.split(":")
                        rng = random.Random(zlib.crc32(("%s|%d|%d|%s|%s" % (setting, rep, k, model, par)).encode()))
                        err = [core.apply_error(model, float(par), rt, tg, spi, rng) for rt, tg in base]
                        est = core.run_apro([core.encode(r, t) for r, t in err], [r for r, _ in err], fixed=True)
                    else:
                        est = method_tree(m, plain, rts, species, td)
                    fn, fp, nt, ne = core.rf_error(est, true_sp)
                    rec.update({"FN": fn, "FP": fp, "nI": nt})
                except Exception as e:  # noqa: BLE001
                    rec.update({"error": "%s: %s" % (type(e).__name__, str(e)[:200])})
                rec["sec"] = round(time.time() - t0, 2)
                with open(out, "a") as f:
                    f.write(json.dumps(rec) + "\n")
                print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
