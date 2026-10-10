"""Core harness for the ASTRAL-Pro rooting/tagging-error experiments.

* simulate GDL families with true root and true D/S tags (gdlsim.py, from the
  sibling claude/cs581-gdl session)
* error models on (rooted gene tree, tags)
* run ASTRAL-Pro3 on given root + tags (patched binary, APRO_FIXED=1, tags encoded
  in branch lengths) or with its own rooting/tagging (unmodified binary)
"""
import os
import random
import subprocess
import tempfile

from gdlsim import simulate
from phylo import Tree, parse_newick, rf_error
import methods as M

APRO_FIXED = os.environ.get("APRO_FIXED_BIN", "/opt/tools/ASTER/bin/astral-pro3-fixed")
APRO = os.environ.get("ASTRALPRO", "/opt/mm/root/envs/gdl/bin/astral-pro3")
ASTRAL = os.environ.get("ASTRAL4", "/opt/mm/root/envs/gdl/bin/astral4")


def sp_of(label):
    return label.rsplit("_", 1)[0]


# ------------------------------------------------------------------ trees + tags

def true_rt(nwk, sp_index):
    t = parse_newick(nwk)
    rt, tags, leafsp = M.root_and_tag(t, sp_of, sp_index, truetags=True)
    return rt, tags


def overlap_tags(rt, sp_index):
    tags, _ = M.tag(rt, sp_of, sp_index)
    return tags


def unrooted_edges(t):
    nb = M._unrooted_adj(t)
    edges = [(u, v) for u in range(len(nb)) for v in nb[u] if u < v]
    return nb, edges


def reroot_random(t, rng):
    nb, edges = unrooted_edges(t)
    a, b = edges[rng.randrange(len(edges))]
    return M._build_rooted(t, nb, a, b)


def clusters(rt):
    """frozenset of leaf labels below each internal node -> node id."""
    below = [None] * len(rt.parent)
    out = {}
    for v in rt.postorder():
        if not rt.children[v]:
            below[v] = frozenset([rt.label[v]])
        else:
            s = frozenset().union(*[below[c] for c in rt.children[v]])
            below[v] = s
            out[s] = v
    return out


def apply_error(model, par, rt, tags, sp_index, rng):
    """Returns (rooted tree, tags). Models:
    true      true root, true tags
    ovl       true root, species-overlap tags (hidden paralogs become S)
    rovl      w.p. p re-root at a uniformly random edge, then overlap tags
    rtrue     w.p. p re-root at a random edge + overlap tags, else true root+true tags
    flip      true root, every internal tag flipped independently w.p. q
    d2s       true root, every true D relabelled S w.p. q
    s2d       true root, every true S relabelled D w.p. q
    """
    if model == "true":
        return rt, tags
    if model == "ovl":
        return rt, overlap_tags(rt, sp_index)
    if model in ("rovl", "rtrue"):
        if rng.random() < par:
            nr = reroot_random(rt, rng)
            return nr, overlap_tags(nr, sp_index)
        return (rt, overlap_tags(rt, sp_index)) if model == "rovl" else (rt, tags)
    if model.startswith("rsp-"):
        # systematic misrooting: w.p. p root on the leaf edge of a random copy of one species
        sp = model[4:]
        cands = [v for v in rt.leaves() if sp_of(rt.label[v]) == sp]
        if cands and rng.random() < par:
            nb = M._unrooted_adj(rt)
            v = rng.choice(cands)
            nr = M._build_rooted(rt, nb, v, nb[v][0])
            return nr, overlap_tags(nr, sp_index)
        return rt, tags
    if model in ("flip", "d2s", "s2d"):
        nt = list(tags)
        for v in range(len(nt)):
            if not rt.children[v]:
                continue
            if model == "flip" or (model == "d2s" and nt[v]) or (model == "s2d" and not nt[v]):
                if rng.random() < par:
                    nt[v] = not nt[v]
        return rt, nt
    raise ValueError(model)


def encode(rt, tags):
    """Newick for the patched binary: internal non-root-child v has length 1 + tag(v);
    leaves 1; root children a, b: len(a) = 1 + tag(a), len(b) = 1 + 2 tag(b) + 8 tag(root)."""
    r = rt.root
    assert len(rt.children[r]) == 2, "root must be binary"

    def rec(u, ln):
        if not rt.children[u]:
            return "%s:%d" % (rt.label[u], ln)
        return "(" + ",".join(rec(c, 1 + int(tags[c] and bool(rt.children[c])))
                              for c in rt.children[u]) + "):%d" % ln
    a, b = rt.children[r]
    la = 1 + int(bool(rt.children[a]) and tags[a])
    lb = 1 + 2 * int(bool(rt.children[b]) and tags[b]) + 8 * int(tags[r])
    return "(" + rec(a, la) + "," + rec(b, lb) + ");"


def plain(rt):
    return rt.newick()


# ------------------------------------------------------------------ runners

def write_mapping(path, trees_rt):
    seen = set()
    with open(path, "w") as f:
        for rt in trees_rt:
            for v in rt.leaves():
                lab = rt.label[v]
                if lab not in seen:
                    seen.add(lab)
                    f.write("%s %s\n" % (lab, sp_of(lab)))


def run_apro(newicks, rts, fixed, threads=1, extra=(), workdir=None):
    with tempfile.TemporaryDirectory(dir=workdir) as td:
        gi, mp, out = (os.path.join(td, x) for x in ("g.nwk", "map.txt", "s.nwk"))
        with open(gi, "w") as f:
            f.write("\n".join(newicks) + "\n")
        write_mapping(mp, rts)
        env = dict(os.environ)
        if fixed:
            env["APRO_FIXED"] = "1"
        else:
            env.pop("APRO_FIXED", None)
        cmd = [APRO_FIXED if fixed else APRO, "-a", mp, "-t", str(threads), "-o", out] + list(extra) + [gi]
        subprocess.run(cmd, check=True, cwd=td, env=env, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        return parse_newick(open(out).read())


def quartet_scores(newicks, rts, fixed, species_nwk, workdir=None):
    """Scores of the three topologies around each branch of species_nwk (freqQuad.csv)."""
    with tempfile.TemporaryDirectory(dir=workdir) as td:
        gi, mp, st, out = (os.path.join(td, x) for x in ("g.nwk", "map.txt", "sp.nwk", "s.nwk"))
        with open(gi, "w") as f:
            f.write("\n".join(newicks) + "\n")
        with open(st, "w") as f:
            f.write(species_nwk + "\n")
        write_mapping(mp, rts)
        env = dict(os.environ)
        if fixed:
            env["APRO_FIXED"] = "1"
        else:
            env.pop("APRO_FIXED", None)
        subprocess.run([APRO_FIXED if fixed else APRO, "-a", mp, "-c", st, "-C", "-u", "3", "-o", out, gi],
                       check=True, cwd=td, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(td, "freqQuad.csv"))]
        return [(r[0], r[1], r[2], float(r[4])) for r in rows]


def fn(est, true_nwk):
    return rf_error(est, parse_newick(true_nwk))[0]
