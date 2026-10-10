"""Species-tree error on 4-taxon pools: fraction of disjoint blocks of K families on which a
method returns a wrong unrooted topology.  Restartable jsonl (method, K, block).
usage: python runmeth.py POOL K NBLOCKS methods(comma) OUT.jsonl"""
import itertools, json, os, re, subprocess, sys, tempfile
from collections import defaultdict
from gdlsim import GNode
from phylo import parse_newick
import apro

AST = "/opt/mm/root/envs/gdl/bin/"
TRUE = "AB|CD"
ROOTED = {}   # all 15 rooted 4-taxon trees
for x, y in [("A", "B"), ("A", "C"), ("A", "D"), ("B", "C"), ("B", "D"), ("C", "D")]:
    rest = [s for s in "ABCD" if s not in (x, y)]
    ROOTED["((%s,%s),%s),%s" % (x, y, rest[0], rest[1])] = (((x, y), rest[0]), rest[1])
    ROOTED["((%s,%s),%s),%s" % (x, y, rest[1], rest[0])] = (((x, y), rest[1]), rest[0])
for (x, y) in [("A", "B"), ("A", "C"), ("A", "D")]:
    rest = tuple(s for s in "ABCD" if s not in (x, y))
    ROOTED["(%s,%s),(%s,%s)" % (x, y, *rest)] = ((x, y), rest)


def topo_of_rooted(t):
    def leaves(z):
        return [z] if isinstance(z, str) else sum((leaves(c) for c in z), [])
    for c in t:
        if not isinstance(c, str) and len(leaves(c)) == 2:
            p = set(leaves(c)); break
    else:
        p = set(leaves([c for c in t if not isinstance(c, str)][0])) - set()
        p = set(leaves(t[0][0])) if not isinstance(t[0][0], str) else None
    side = p if "A" in p else set("ABCD") - p
    return "".join(sorted(side)) + "|" + "".join(sorted(set("ABCD") - side))


def topo_of_newick(nwk):
    t = parse_newick(nwk)
    below = {}
    for v in t.postorder():
        below[v] = {t.label[v].split("_")[0]} if not t.children[v] else set().union(*[below[c] for c in t.children[v]])
        if len(below[v]) == 2:
            side = below[v] if "A" in below[v] else set("ABCD") - below[v]
            return "".join(sorted(side)) + "|" + "".join(sorted(set("ABCD") - side))
    raise ValueError(nwk)


def parse_g(nwk):
    """newick with D/S internal labels -> GNode"""
    t = parse_newick(nwk)
    def rec(v):
        if not t.children[v]:
            return GNode("L", label=t.label[v])
        return GNode(t.label[v] or "S", [rec(c) for c in t.children[v]])
    return rec(t.root)


def plain(nwk):
    return re.sub(r"\)[DS]", ")", nwk)


def argmax(sc):
    best = max(sc.values()) if sc else 0
    win = sorted(k for k, v in sc.items() if v == best)
    return win[0] if len(win) == 1 else "tie"


def py_method(m, gs, own_cache):
    tot = defaultdict(float)
    for i, g in enumerate(gs):
        if m == "true":
            r = apro.truetag_scores(g)
        elif m == "ovl":
            r = apro.family_scores(g, "true")
        elif m == "own":
            r = apro.family_scores(g, "own")
        elif m == "recon_true":
            r = apro.recon_scores(g, (((("A", "B"), "C"), "D")))
        elif m == "recon_first":
            r = apro.recon_scores(g, own_cache["S0"])
        for k, v in r.items():
            if "|" in k:
                tot[k] += v
    return argmax(tot), dict(tot)


def gtp(gs, loss_w):
    cost = {}
    for name, st in ROOTED.items():
        cost[name] = sum(apro.gtp_cost(g, st, loss_w) for g in gs)
    b = min(cost.values())
    win = [n for n, c in cost.items() if c == b]
    tops = {topo_of_rooted(ROOTED[n]) for n in win}
    return (tops.pop() if len(tops) == 1 else "tie"), {}


def binary_method(m, nwks, td):
    gi = os.path.join(td, "g.nwk"); mp = os.path.join(td, "map.txt"); out = os.path.join(td, "s.nwk")
    pl = [plain(x) for x in nwks]
    open(gi, "w").write("\n".join(pl) + "\n")
    labs = sorted(set(re.findall(r"[A-D]_\d+", "".join(pl))))
    open(mp, "w").write("".join("%s %s\n" % (l, l.split("_")[0]) for l in labs))
    q = dict(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True, cwd=td)
    if m == "apro_bin":
        subprocess.run([AST + "astral-pro3", "-a", mp, "-o", out, gi], **q)
    elif m == "astral_multi":
        subprocess.run([AST + "astral4", "-a", mp, "-o", out, gi], **q)
    elif m == "disco_astral":
        dd = os.path.join(td, "disco.nwk")
        subprocess.run(["python3", "/opt/src/DISCO/disco.py", "-i", gi, "-o", dd, "-d", "_"], **q)
        dl = [l for l in open(dd) if ";" in l]
        labs2 = sorted(set(re.findall(r"[A-D]_\d+", "".join(dl))))
        open(mp, "w").write("".join("%s %s\n" % (l, l.split("_")[0]) for l in labs2))
        subprocess.run([AST + "astral4", "-a", mp, "-o", out, dd], **q)
    elif m == "wqfm_gdl":
        subprocess.run(["bash", "/opt/src/wQFM-GDL/wQFM-GDL.sh", "-i", gi, "-o", out, "-q"], **q)
    elif m == "duploss2":
        sp = [re.sub(r"([A-D])_\d+", r"\1", x) for x in pl]
        gi2 = os.path.join(td, "g2.nwk"); open(gi2, "w").write("\n".join(sp) + "\n")
        subprocess.run(["/opt/src/DupLoss-2/Executables/DupLoss-2.linux", "-i", gi2, "-o", out], **q)
    txt = open(out).read()
    nw = [l for l in txt.splitlines() if ";" in l][-1]
    nw = nw[nw.index("("):]
    return topo_of_newick(nw), {}


if __name__ == "__main__":
    pool, K, NB, methods, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4].split(","), sys.argv[5]
    lines = [l.strip() for l in open(pool) if ";" in l]
    done = set()
    if os.path.exists(out):
        done = {(r["method"], r["K"], r["block"]) for r in map(json.loads, open(out))}
    for b in range(min(NB, len(lines) // K)):
        blk = lines[b * K:(b + 1) * K]
        gs = None
        cache = {}
        for m in methods:
            if (m, K, b) in done:
                continue
            if gs is None:
                gs = [parse_g(x) for x in blk]
            if m in ("true", "ovl", "own", "recon_true", "recon_first"):
                if m == "recon_first":
                    s0, _ = py_method("own", gs, cache)
                    side = s0.split("|")[0] if s0 != "tie" else "AB"
                    other = [c for c in "ABCD" if c not in side]
                    cache["S0"] = (((side[0], side[1]), other[0]), other[1])
                res, sc = py_method(m, gs, cache)
            elif m in ("gtp_dup", "gtp_dl"):
                res, sc = gtp(gs, 0.0 if m == "gtp_dup" else 1.0)
            else:
                with tempfile.TemporaryDirectory() as td:
                    try:
                        res, sc = binary_method(m, blk, td)
                    except Exception as e:
                        res, sc = "ERR:" + str(e)[:80], {}
            with open(out, "a") as f:
                f.write(json.dumps({"pool": os.path.basename(pool), "method": m, "K": K, "block": b, "topo": res,
                                    "wrong": res != TRUE, "scores": sc}) + "\n")
            print(b, m, res, flush=True)
