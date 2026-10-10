"""Pure-GDL simulator that keeps branch lengths (time units), for sequence simulation.
Species tree: phylo.Tree with lengths; lam/mu dicts by node id (rate on branch above).
Family starts with one copy at the species root (no root branch). Leaves 'Sp_k'."""
import random
from phylo import parse_newick


def sim_family(st, lam, mu, rng, cap=3000, tags=False):
    cnt, box = {}, [cap]

    def along(v, rem, acc):
        """copy entering branch above v with remaining time rem; acc = length so far."""
        l, m = lam[v], mu[v]
        tot = l + m
        while True:
            dt = rng.expovariate(tot) if tot > 0 else float("inf")
            if dt >= rem:
                return at_node(v, acc + rem)
            rem -= dt; acc += dt
            if rng.random() < l / tot:
                box[0] -= 1
                if box[0] < 0:
                    raise OverflowError
                a, b = along(v, rem, 0.0), along(v, rem, 0.0)
                return join(a, b, acc, "D")
            return None

    def join(a, b, acc, kind="S"):
        if a is None and b is None:
            return None
        if a is None:
            return (b[0], b[1] + acc)
        if b is None:
            return (a[0], a[1] + acc)
        return ("(%s:%.6f,%s:%.6f)%s" % (a[0], a[1], b[0], b[1], kind if tags else ""), acc)

    def at_node(v, acc):
        ch = st.children[v]
        if not ch:
            sp = st.label[v]; k = cnt.get(sp, 0); cnt[sp] = k + 1
            return ("%s_%d" % (sp, k), acc)
        kids = [along(c, st.length[c], 0.0) for c in ch]
        kids = [k for k in kids if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return (kids[0][0], kids[0][1] + acc)
        return join(kids[0], kids[1], acc)

    try:
        r = at_node(st.root, 0.0)
    except OverflowError:
        return None
    if r is None or "," not in r[0]:
        return None
    return r[0] + ";", cnt
