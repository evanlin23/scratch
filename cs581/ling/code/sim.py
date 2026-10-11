"""Simulator of linguistic characters on a (dated) tree with optional contact (borrowing) edges.

Model (in the spirit of Warnow-Evans-Ringe-Nakhleh 2006, plus polymorphism):
* Model tree: Yule topology, root at time 0, contemporaneous leaves at time 1; a fraction
  `anc` of leaves are 'ancient' (terminal edge cut short). Edge rate multipliers are
  lognormal(0, dlc) (deviation from the lexical clock).
* Character rates: type base rate x Gamma(alpha) across characters x lognormal(0, het) per
  (character, edge) (heterotachy / non-i.i.d. evolution).
* Lexical characters (meaning slots): every change creates a NEW cognate class (no back
  mutation) except in the fraction `hom_L` of homoplastic characters, where a change goes
  w.p. 1/2 to a state already used (parallel semantic shift). Polymorphism: a change on a
  monomorphic slot adds a synonym w.p. `p_poly`; on a polymorphic slot it is a loss of one
  word w.p. `q_loss`, else a replacement of one word.
* Morphological / phonological characters: fraction `hom_M` / `hom_P` homoplastic. Homoplastic
  phonological characters are binary (0<->1: parallel and back mutation, i.e. simple natural
  sound changes); homoplastic morphological characters revert to a previously used state w.p.
  1/2. Others are infinite-state. No polymorphism.
* Borrowing: `n_contact` contact edges between lineages that coexist at a random time; over a
  window of length `contact_len` each character of type t is borrowed (donor -> recipient,
  fixed direction per contact edge) w.p. `b_<t>`. A borrowed lexical item replaces the
  recipient's word, or w.p. `p_poly` is added as a synonym.

Output: list of per-character columns of frozensets (leaf state sets), types, and the true tree.
"""
import numpy as np

from trees import ArrTree

DEFAULTS = dict(n=24, nL=250, nM=15, nP=25,
                rL=4.0, rM=2.0, rP=0.6,          # base rates (expected changes per unit height)
                alpha=1.0, dlc=0.3, het=0.3, anc=0.3,
                hom_L=0.10, hom_M=0.20, hom_P=0.30,
                p_poly=0.3, q_loss=0.6,
                n_contact=1, contact_len=0.15, b_L=0.10, b_M=0.0, b_P=0.03)


def yule_tree(n, rng):
    """Returns ArrTree with node times (root=0, contemporaneous leaves=1)."""
    left = -np.ones(2 * n - 1, np.int64); right = left.copy(); parent = left.copy()
    time = np.zeros(2 * n - 1)
    # grow lineages forward in time; tips are open lineages identified by temporary ids
    nxt_int = 2 * n - 2
    root = nxt_int; nxt_int -= 1
    tips = [(root, 0), (root, 1)]  # (parent internal node, side)
    t = 0.0
    while len(tips) < n:
        t += rng.exponential(1.0 / len(tips))
        i = rng.integers(len(tips))
        par, side = tips.pop(i)
        v = nxt_int; nxt_int -= 1
        time[v] = t
        parent[v] = par
        if side == 0:
            left[par] = v
        else:
            right[par] = v
        tips += [(v, 0), (v, 1)]
    T = t + rng.exponential(1.0 / len(tips))
    order = rng.permutation(n)
    for (par, side), leaf in zip(tips, order):
        parent[leaf] = par
        if side == 0:
            left[par] = leaf
        else:
            right[par] = leaf
        time[leaf] = T
    assert nxt_int == n - 1
    time /= T
    t = ArrTree(n, left, right, parent, root)
    t.time = time
    return t


def simulate(params=None, rng=None, tree=None):
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    rng = rng or np.random.default_rng()
    n = p['n']
    t = tree if tree is not None else yule_tree(n, rng)
    N = 2 * n - 1
    time = t.time.copy()
    # ancient languages: shorten terminal edges
    for v in range(n):
        if rng.random() < p['anc']:
            el = time[v] - time[t.parent[v]]
            time[v] -= rng.uniform(0, 0.6) * el
    edge_mult = rng.lognormal(0, p['dlc'], N)
    edge_mult /= edge_mult.mean()

    types = ['L'] * p['nL'] + ['M'] * p['nM'] + ['P'] * p['nP']
    C = len(types)
    tarr = np.array(types)
    base = np.where(tarr == 'L', p['rL'], np.where(tarr == 'M', p['rM'], p['rP']))
    crate = base * rng.gamma(p['alpha'], 1.0 / p['alpha'], C)
    hetm = rng.lognormal(0, p['het'], (C, N))
    hom = np.zeros(C, bool)
    for ty in 'LMP':
        idx = np.where(tarr == ty)[0]
        k = int(round(p['hom_' + ty] * len(idx)))
        hom[rng.choice(idx, k, replace=False)] = True
    binary = hom & (tarr == 'P')
    is_lex = tarr == 'L'

    s1 = {}; s2 = {}; cur = {}
    nxt_state = np.ones(C, np.int64)
    pools = [[0] for _ in range(C)]
    stats = dict(changes=0, borrow=0, poly_events=0)

    def new_state(c):
        x = int(nxt_state[c]); nxt_state[c] += 1
        pools[c].append(x)
        return x

    def mutate(e, c):
        a = s1[e]; b = s2[e]
        stats['changes'] += 1
        if binary[c]:
            a[c] = 1 - a[c]
            return
        if hom[c] and len(pools[c]) > 1 and rng.random() < 0.5:
            cand = [x for x in pools[c] if x != a[c] and x != b[c]]
            newv = cand[rng.integers(len(cand))] if cand else new_state(c)
        else:
            newv = new_state(c)
        if is_lex[c]:
            if b[c] < 0:
                if rng.random() < p['p_poly']:
                    b[c] = newv; stats['poly_events'] += 1
                else:
                    a[c] = newv
            else:
                if rng.random() < p['q_loss']:
                    if rng.random() < 0.5:
                        a[c] = b[c]
                    b[c] = -1
                else:
                    if rng.random() < 0.5:
                        a[c] = newv
                    else:
                        b[c] = newv
        else:
            a[c] = newv

    def advance(e, to):
        dt = to - cur[e]
        if dt <= 0:
            return
        lam = crate * edge_mult[e] * hetm[:, e] * dt
        k = rng.poisson(lam)
        for c in np.nonzero(k)[0]:
            for _ in range(k[c]):
                mutate(e, c)
        cur[e] = to

    # contact edges -> borrowing events
    events = []  # (time, kind, payload)
    contacts = []
    for _ in range(p['n_contact']):
        for _try in range(100):
            tc = rng.uniform(0.25, 0.9)
            alive = [v for v in range(N) if v != t.root and time[t.parent[v]] < tc < time[v]]
            if len(alive) < 2:
                continue
            d, r = rng.choice(alive, 2, replace=False)
            if t.parent[d] == t.parent[r]:
                continue  # sisters: borrowing undetectable-ish, skip
            end = min(tc + p['contact_len'], time[d], time[r])
            contacts.append((int(d), int(r), tc, end))
            for c in range(C):
                if rng.random() < p['b_' + types[c]]:
                    events.append((rng.uniform(tc, end), 1, (int(d), int(r), c)))
            break
    for v in range(n, N):
        events.append((time[v], 0, v))
    events.sort(key=lambda x: (x[0], x[1]))

    root = t.root
    s1[root] = np.zeros(C, np.int64); s2[root] = -np.ones(C, np.int64); cur[root] = 0.0
    for tm, kind, pay in events:
        if kind == 0:
            v = pay
            if v != root:
                advance(v, tm)
            for ch in (t.left[v], t.right[v]):
                s1[ch] = s1[v].copy(); s2[ch] = s2[v].copy(); cur[ch] = tm
        else:
            d, r, c = pay
            if d not in cur or r not in cur:
                continue
            advance(d, tm); advance(r, tm)
            word = s1[d][c] if (s2[d][c] < 0 or rng.random() < 0.5) else s2[d][c]
            stats['borrow'] += 1
            if is_lex[c] and s2[r][c] < 0 and rng.random() < p['p_poly'] and s1[r][c] != word:
                s2[r][c] = word
            else:
                s1[r][c] = word; s2[r][c] = -1
    for v in range(n):
        advance(v, time[v])
    chars = []
    for c in range(C):
        col = []
        for v in range(n):
            col.append(frozenset([s1[v][c]]) if s2[v][c] < 0 else frozenset([s1[v][c], s2[v][c]]))
        chars.append(col)
    info = dict(contacts=contacts, hom=hom, stats=stats)
    return t, chars, types, info


def resolve_random(chars, rng):
    """Lexicographer-style coding: polymorphic leaf -> one of its states at random."""
    out = []
    for col in chars:
        out.append([st if len(st) == 1 else frozenset([sorted(st)[rng.integers(len(st))]]) for st in col])
    return out


def drop_polymorphic(chars, types):
    keep = [i for i, col in enumerate(chars) if all(len(s) == 1 for s in col)]
    return [chars[i] for i in keep], [types[i] for i in keep]


def summary(chars, types):
    """Summary stats comparable with the IE data: mean #states per type, % polymorphic."""
    out = {}
    for ty in 'LMP':
        cols = [c for c, t in zip(chars, types) if t == ty]
        if cols:
            out['states_' + ty] = float(np.mean([len(set().union(*c)) for c in cols]))
    lex = [c for c, t in zip(chars, types) if t == 'L']
    out['poly_char_frac'] = float(np.mean([any(len(s) > 1 for s in c) for c in lex])) if lex else 0.0
    out['poly_cell_frac'] = float(np.mean([len(s) > 1 for c in lex for s in c])) if lex else 0.0
    return out
