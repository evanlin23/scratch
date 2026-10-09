"""Minimal Newick parsing and centroid-edge bipartition (no dependencies)."""

import re


def leaf_sets(newick):
    tokens = re.findall(r"\(|\)|,|;|[^(),;]+", newick.strip())
    children, labels = [[]], [None]
    stack, cur, prev = [], 0, None
    for tok in tokens:
        if tok == "(":
            children.append([]); labels.append(None)
            new = len(children) - 1
            children[cur].append(new)
            stack.append(cur)
            cur = new
        elif tok == ")":
            cur = stack.pop()
        elif tok in (",", ";"):
            pass
        elif prev != ")":
            name = tok.split(":")[0].strip().strip("'")
            if name:
                children.append([]); labels.append(name)
                children[cur].append(len(children) - 1)
        prev = tok
    # clade (leaf set) of every node, post-order
    clades = [None] * len(children)
    order, todo = [], [0]
    while todo:
        v = todo.pop()
        order.append(v)
        todo.extend(children[v])
    for v in reversed(order):
        clades[v] = {labels[v]} if labels[v] is not None else set().union(*(clades[c] for c in children[v]))
    return clades


def centroid_split(newick):
    """The edge whose removal gives the most balanced bipartition: (A, B) with |A| <= |B|."""
    clades = leaf_sets(newick)
    every = clades[0]
    n = len(every)
    best = min((c for c in clades[1:] if 0 < len(c) < n), key=lambda c: abs(n - 2 * len(c)))
    a, b = sorted(best), sorted(every - best)
    return (a, b) if len(a) <= len(b) else (b, a)
