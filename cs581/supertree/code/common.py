"""Shared tree utilities for the supertree pilot (treeswift-based, bitmask splits)."""
import os, re, subprocess, time, resource
import treeswift

ENV = "/opt/mamba/envs/st/bin"
ASTRAL3 = "/opt/runs/ASTRAL/Astral/astral.5.7.8.jar"
GTM = "/opt/runs/GTM/gtm.py"


def read_trees(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(treeswift.read_tree_newick(line))
    return out


def leaves(t):
    return [n.label for n in t.traverse_leaves()]


def strip_newick(t):
    """Topology-only newick (no lengths, no internal labels)."""
    for n in t.traverse_preorder():
        n.edge_length = None
        if not n.is_leaf():
            n.label = None
    return t.newick()


def restrict(t, keep):
    """Induced subtree on leaf set `keep` (set of labels); returns a copy or None if <3 leaves."""
    ls = [l for l in leaves(t) if l in keep]
    if len(ls) < 3:
        return None
    r = t.extract_tree_with(set(ls), suppress_unifurcations=True)
    return r


def splits(t, index):
    """Non-trivial bipartitions as canonical int bitmasks over `index` (label->bit).
    Leaves not in index are ignored (caller should restrict first)."""
    full = 0
    for l in leaves(t):
        full |= 1 << index[l]
    nleaf = bin(full).count("1")
    out = set()
    mask = {}
    for n in t.traverse_postorder():
        if n.is_leaf():
            mask[n] = 1 << index[n.label]
        else:
            m = 0
            for c in n.children:
                m |= mask[c]
            mask[n] = m
            k = bin(m).count("1")
            if 2 <= k <= nleaf - 2:
                if m & 1:  # canonical: side not containing bit 0
                    m = full ^ m
                out.add(m)
    return out, full


def fn_fp(true_t, est_t):
    """FN and FP rates (unrooted) on the common leaf set; also RF rate = (FN+FP)/2 counts / (n-3)."""
    common = sorted(set(leaves(true_t)) & set(leaves(est_t)))
    keep = set(common)
    idx = {l: i for i, l in enumerate(common)}
    tt = restrict(true_t, keep) if len(keep) < len(leaves(true_t)) else true_t
    et = restrict(est_t, keep) if len(keep) < len(leaves(est_t)) else est_t
    s_true, _ = splits(tt, idx)
    s_est, _ = splits(et, idx)
    fn = len(s_true - s_est)
    fp = len(s_est - s_true)
    n = len(common)
    return {"n": n, "fn": fn / max(1, len(s_true)), "fp": fp / max(1, len(s_est)) if s_est else 0.0,
            "rf": (fn + fp) / (2 * (n - 3)), "fn_cnt": fn, "fp_cnt": fp,
            "n_true_splits": len(s_true), "n_est_splits": len(s_est)}


def run(cmd, log=None, timeout=None):
    """Run a command, return (wall seconds, peak child RSS MB)."""
    t0 = time.time()
    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    with open(log or os.devnull, "a") as lf:
        p = subprocess.run(cmd, shell=True, stdout=lf, stderr=lf, timeout=timeout)
    wall = time.time() - t0
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if p.returncode != 0:
        raise RuntimeError(f"command failed ({p.returncode}): {cmd}")
    return wall, peak / 1024.0


def timed(cmd, log=None, timeout=None):
    """Run via /usr/bin/time -v to get the peak RSS of this specific command."""
    tf = (log or "/tmp/x") + ".time"
    full = f"/usr/bin/time -v -o {tf} bash -c {sh_quote(cmd)}"
    t0 = time.time()
    with open(log or os.devnull, "a") as lf:
        p = subprocess.run(full, shell=True, stdout=lf, stderr=lf, timeout=timeout)
    wall = time.time() - t0
    rss = None
    try:
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(tf).read())
        rss = int(m.group(1)) / 1024.0 if m else None
    except OSError:
        pass
    if p.returncode != 0:
        raise RuntimeError(f"command failed ({p.returncode}): {cmd}")
    return wall, rss


def sh_quote(s):
    return "'" + s.replace("'", "'\"'\"'") + "'"
