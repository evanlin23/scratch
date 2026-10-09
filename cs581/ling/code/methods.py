"""Tree-estimation methods for multistate linguistic characters.

Input everywhere: `chars` = list over characters of length-n lists of state sets (frozensets of
ints; a set with >1 element is a polymorphic leaf), `types` = list of 'L','M','P' per character.
Each method returns an ArrTree with leaves in the same order as the character columns.
"""
import os
import subprocess
import tempfile
import time

import numpy as np

import fitch
from trees import nj, parse_newick, from_nested

ALPH = '0123456789ABCDEFGHIJKLMNOPQRSTUV'


def hamming(chars, n):
    """Pairwise distance = fraction of characters on which two languages share no state."""
    D = np.zeros((n, n))
    for col in chars:
        for i in range(n):
            for j in range(i + 1, n):
                if not (col[i] & col[j]):
                    D[i, j] += 1
    D = D + D.T
    return D / max(1, len(chars))


def nj_tree(chars, n, correct=True):
    p = hamming(chars, n)
    if correct:
        # infinite-alleles style correction d = -ln(1-p), capped
        D = -np.log(np.clip(1 - p, 0.02, 1))
    else:
        D = p
    return nj(D)


def char_search(chars, types, n, kind, weights=None, cap=None, rng=None, nstarts=3, start_trees=()):
    ls, lb = fitch.encode_sets(chars, n)
    w = np.ones(len(chars)) if weights is None else np.array([weights[t] for t in types], float)
    obj = fitch.Objective(kind, cap=cap)
    starts = list(start_trees) + ['random'] * nstarts
    t, sc = fitch.search(ls, lb, w, obj, starts, rng)
    return t, sc


def _run_iqtree(phy_text, args, workdir, seed):
    path = os.path.join(workdir, 'aln.phy')
    with open(path, 'w') as f:
        f.write(phy_text)
    cmd = ['iqtree2', '-s', path, '-seed', str(seed), '-nt', '1', '--quiet', '-redo'] + args
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(path + '.treefile') as f:
        return f.read()


def _names(n):
    return ['t%d' % i for i in range(n)]


def _nested_to_tree(nwk, n):
    names = _names(n)
    nested = parse_newick(nwk)
    return from_nested(nested, names)


def iqtree_mk(chars, n, seed=1, poly='random', rng=None):
    """Multistate Mk+G4 (IQ-TREE MORPH). Polymorphic leaves -> one state at random."""
    rng = rng or np.random.default_rng(seed)
    cols = []
    for col in chars:
        lab = {}
        out = []
        for st in col:
            x = sorted(st)[rng.integers(len(st))] if len(st) > 1 else next(iter(st))
            if x not in lab:
                lab[x] = len(lab)
            out.append(lab[x])
        if len(lab) > len(ALPH):
            continue
        cols.append(out)
    names = _names(n)
    lines = ['%d %d' % (n, len(cols))]
    for i in range(n):
        lines.append(names[i] + ' ' + ''.join(ALPH[c[i]] for c in cols))
    with tempfile.TemporaryDirectory() as d:
        nwk = _run_iqtree('\n'.join(lines) + '\n', ['-st', 'MORPH', '-m', 'MK+G4'], d, seed)
    return _nested_to_tree(nwk, n)


def binary_matrix(chars, n):
    """Gray & Atkinson style binary encoding: one presence/absence column per state."""
    cols = []
    for col in chars:
        states = sorted(set().union(*col))
        if len(states) < 2:
            continue
        for s in states:
            cols.append([1 if s in st else 0 for st in col])
    return cols


def iqtree_binary(chars, n, seed=1):
    cols = binary_matrix(chars, n)
    names = _names(n)
    lines = ['%d %d' % (n, len(cols))]
    for i in range(n):
        lines.append(names[i] + ' ' + ''.join(str(c[i]) for c in cols))
    with tempfile.TemporaryDirectory() as d:
        nwk = _run_iqtree('\n'.join(lines) + '\n', ['-st', 'BIN', '-m', 'GTR2+FO+G4'], d, seed)
    return _nested_to_tree(nwk, n)


def timed(fn, *a, **k):
    t0 = time.time()
    out = fn(*a, **k)
    return out, time.time() - t0
