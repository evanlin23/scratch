"""Figures of the counterexamples for theory.md, each verified numerically (equal d vectors).
  fig_claw.png      : k = 3 closed claw, n = 5, generic lengths: three different exact solutions
  fig_adjacent.png  : k = 3 with an adjacent pair (the paper's Fig. 1d situation) with GENERIC lengths
  fig_splitclaw.png : k = 4 split claw, n = 7 (non-adjacent, no closed claw, still a continuum)
"""
import sys, os, itertools, collections
import numpy as np
import networkx as nx
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "decodiphy", "code"))
from config_enum import laplacian
OUT = os.path.join(os.path.dirname(__file__), "..", "figs")
COL = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]


def dist(N, E, lens):
    G = nx.Graph()
    for (a, b), l in zip(E, lens): G.add_edge(a, b, weight=l)
    D = dict(nx.all_pairs_dijkstra_path_length(G))
    return np.array([[D[u][v] for v in range(N)] for u in range(N)])


def dvec(n, E, lens, sol):
    """sol: list of (edge (a,b), p, x from a), ybar"""
    pls, y = sol
    N = 1 + max(max(e) for e in E)
    D = dist(N, E, lens)
    L = dict(zip(E, lens)); L.update({(b, a): l for (a, b), l in zip(E, lens)})
    m = np.zeros(N)
    for (a, b), p, x in pls:
        m[a] += p * (1 - x); m[b] += p * x
    return D[:n] @ m + y, m


def measure_to_solution(E, S, m):
    """realize node measure m on edge set S: split shared nodes equally, then p, x per edge."""
    cnt = collections.Counter(u for e in S for u in e)
    out = []
    for a, b in S:
        sa, sb = m[a] / cnt[a], m[b] / cnt[b]
        out.append(((a, b), sa + sb, sb / (sa + sb)))
    return out


def draw(ax, n, E, lens, pos, sol, title, col):
    pls, y = sol
    for (a, b) in E:
        ax.plot(*zip(pos[a], pos[b]), color="0.6", lw=2, zorder=1)
    for u in pos:
        if u < n:
            ax.text(*pos[u], f" {u}", fontsize=9, va="center")
        ax.plot(*pos[u], "o", color="k" if u < n else "0.4", ms=4, zorder=2)
    for j, ((a, b), p, x) in enumerate(pls):
        P = (1 - x) * np.array(pos[a]) + x * np.array(pos[b])
        ax.plot(*P, "s", color=col, ms=6 + 18 * p, alpha=.85, zorder=3)
        ax.text(P[0], P[1] + .12, f"p={p:.2f}\nx={x:.2f}", fontsize=7, ha="center", color=col)
    ax.set_title(title + f"\n$\\bar y$={y:.3f}", fontsize=9); ax.axis("off"); ax.set_aspect("equal")


def fig_claw():
    n = 5
    E = [(2, 5), (5, 6), (6, 0), (6, 3), (5, 7), (7, 1), (7, 4)]
    lens = [0.7, 0.4, 0.9, 0.5, 0.3, 0.8, 0.6]
    pos = {5: (0, 0), 2: (0, -1.3), 6: (-1, .6), 0: (-2, 1.4), 3: (-2, -.2), 7: (1, .6), 1: (2, 1.4), 4: (2, -.2)}
    sols = [([((2, 5), .2, .8), ((6, 0), .4, .5), ((7, 1), .4, .5)], .1)]
    d0, m0 = dvec(n, E, lens, sols[0])
    L = laplacian(8, E, lens)
    # moves alpha = eps * delta_6 (claw centred at 6: 6 and its neighbours 5,0,3 are in V(S))
    for eps, S in ((-0.0281, [(5, 6), (6, 0), (7, 1)]), (None, None)):
        pass
    # exact alternatives reported by the exhaustive LP search of the first pilot
    sols += [([((5, 6), .4283, .1518), ((6, 0), .265, .7547), ((7, 1), .3067, .6522)], .128),
             ([((6, 0), .33, .6061), ((5, 7), .4167, .128), ((7, 1), .2533, .7895)], .128)]
    errs = [np.abs(dvec(n, E, lens, s)[0] - d0).max() for s in sols]
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))
    for ax, s, c, t in zip(axs, sols, COL, ["truth", "alternative 1 (edge (2,5) -> (5,6))", "alternative 2 (edge (2,5) -> (5,7))"]):
        draw(ax, n, E, lens, pos, s, t, c)
    fig.suptitle(f"k = 3, closed claws at 6 and 7 (generic lengths); max |d - d_truth| = {max(errs):.1e} (rounding of reported p,x)", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_claw.png", dpi=130)
    return errs


def alternative_from_direction(n, E, lens, S, Sp, alpha_int, rng):
    """build m on V(S), m' on V(S') with m' - m = -L alpha (scaled)."""
    N = 2 * n - 2
    L = laplacian(N, E, lens)
    g = -L[:, n:] @ alpha_int
    VS = {u for e in S for u in e}; VSp = {u for e in Sp for u in e}
    mp = np.zeros(N)
    for u in VSp:
        mp[u] = rng.uniform(1, 2) + (abs(g[u]) if u in VS else 0)
    s = 1.0
    for u in VSp - VS:
        s = min(s, mp[u] / g[u])  # need m_u = mp_u - s g_u = 0 -> choose s exactly per node below
    # choose s so that nodes in V(S')\V(S) vanish: requires mp_u = s g_u; set mp accordingly
    s = 1.0
    for u in VSp - VS:
        mp[u] = s * g[u]
    m = mp - s * g
    assert np.all(m[list(VS)] > 0) and np.all(np.abs(m[[u for u in range(N) if u not in VS]]) < 1e-12), m
    tot = m.sum(); m /= tot; mp /= tot
    ybar = 0.2
    ybarp = ybar - s * alpha_int.sum() / tot
    return m, mp, ybar, ybarp


def lp_alpha(L, n, VS, VSp):
    from scipy.optimize import linprog
    N = L.shape[0]; G = -L[:, n:]
    zero = [u for u in range(N) if u not in VS and u not in VSp]
    pos = [u for u in VSp if u not in VS]; neg = [u for u in VS if u not in VSp]
    if not pos and not neg:
        return None
    A_ub = np.vstack([-G[pos], G[neg]]); b_ub = -np.ones(len(pos) + len(neg))
    r = linprog(np.zeros(N - n), A_ub=A_ub, b_ub=b_ub, A_eq=G[zero] if zero else None,
                b_eq=np.zeros(len(zero)) if zero else None, bounds=[(-10, 10)] * (N - n), method="highs")
    return r.x if r.status == 0 else None


def fig_adjacent():
    # paper Fig 1d shape: caterpillar ((1,2)A, 3 B, (4,5)C) -> leaves 0..4 = 1..5; internal 5=A,6=B,7=C
    n = 5
    E = [(0, 5), (1, 5), (5, 6), (2, 6), (6, 7), (3, 7), (4, 7)]
    rng = np.random.default_rng(11)
    lens = [float(v) for v in np.round(rng.uniform(0.3, 1.0, 7), 2)]
    pos = {0: (-2, 1), 1: (-2, -1), 5: (-1, 0), 6: (0, 0), 2: (0, -1.3), 7: (1, 0), 3: (2, 1), 4: (2, -1)}
    S = [(0, 5), (1, 5), (2, 6)]
    L = laplacian(8, E, lens)
    VS = {u for e in S for u in e}
    found = []
    for Sp in itertools.combinations(E, 3):
        if set(Sp) == set(S):
            continue
        VSp = {u for e in Sp for u in e}
        a = lp_alpha(L, n, VS, VSp)
        if a is not None:
            found.append((list(Sp), a))
    Sp, alpha = found[0]
    g = -L[:, n:] @ alpha
    g /= np.abs(g).max(); alpha = alpha / np.abs(-L[:, n:] @ alpha).max()
    VSp = {u for e in Sp for u in e}
    m = np.zeros(8)
    for u in VS:
        m[u] = rng.uniform(0.5, 1.0)
    s = 0.3
    for u in VS - VSp:
        m[u] = -s * g[u]          # these vanish in m' = m + s g
    mp = m + s * g
    assert np.all(mp[list(VSp)] > 0) and np.all(np.abs(mp[[u for u in range(8) if u not in VSp]]) < 1e-12)
    tot = m.sum(); m /= tot; mp /= tot
    y = 0.2; yp = y - s * alpha.sum() / tot
    sol = (measure_to_solution(E, S, m), y)
    solp = (measure_to_solution(E, Sp, mp), yp)
    d0, _ = dvec(n, E, lens, sol); d1, _ = dvec(n, E, lens, solp)
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.4))
    draw(axs[0], n, E, lens, pos, sol, "truth: adjacent pair at node 5", COL[0])
    draw(axs[1], n, E, lens, pos, solp, f"alternative S'={Sp}", COL[1])
    fig.suptitle(f"paper's Fig. 1d shape with random lengths {lens}; max |d-d'| = {np.abs(d0-d1).max():.1e}", fontsize=8)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_adjacent.png", dpi=130)
    return dict(err=float(np.abs(d0 - d1).max()), S=S, Sp=Sp, n_alternatives=len(found), lens=lens, sol=sol, solp=solp)


def fig_splitclaw():
    n = 7  # leaves 0..6; internal 7=v, 8=u, 9=w, 10, 11
    E = [(7, 0), (7, 10), (7, 8), (8, 6), (8, 9), (9, 2), (9, 11), (10, 1), (10, 3), (11, 4), (11, 5)]
    rng = np.random.default_rng(5)
    lens = [float(v) for v in np.round(rng.uniform(0.3, 1.0, len(E)), 2)]
    pos = {7: (-1, 0), 8: (0, 0), 9: (1, 0), 6: (0, -1.2), 0: (-1.6, -1), 10: (-1.8, 0.8), 1: (-2.8, 0.4), 3: (-2.2, 1.8),
           2: (1.6, -1), 11: (1.8, 0.8), 4: (2.8, 0.4), 5: (2.2, 1.8)}
    S = [(7, 0), (10, 1), (9, 2), (11, 4)]
    N = 12
    L = laplacian(N, E, lens)
    VS = {u for e in S for u in e}
    out = [u for u in range(N) if u not in VS]
    _, sv, Vt = np.linalg.svd(L[out][:, n:])
    alpha = Vt[-1]; assert np.abs(L[out][:, n:] @ alpha).max() < 1e-10
    g = -L[:, n:] @ alpha
    assert np.all(np.abs(g[out]) < 1e-10)
    m = np.zeros(N)
    for u in VS:
        m[u] = rng.uniform(0.5, 1.0)
    m /= m.sum()
    # feasible range of t: m + t g >= 0 on V(S)
    ts = [-m[u] / g[u] for u in VS if abs(g[u]) > 1e-12]
    tlo = max(t for t in ts if t < 0); thi = min(t for t in ts if t > 0)
    y = 0.3
    sols = []
    for t in (0, 0.9 * tlo, 0.9 * thi):
        mt = m + t * g
        sols.append((measure_to_solution(E, S, mt), y - t * alpha.sum()))
    d = [dvec(n, E, lens, s)[0] for s in sols]
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.8))
    for ax, s, c, ttl in zip(axs, sols, COL, ["truth", "same edges, t<0", "same edges, t>0"]):
        draw(ax, n, E, lens, pos, s, ttl, c)
    err = max(np.abs(d[1] - d[0]).max(), np.abs(d[2] - d[0]).max())
    fig.suptitle(f"k = 4 split claw: v=7, w=9 in V(S), u=8 not; alpha supported on {{7, 9}}; max |d - d_truth| = {err:.1e}", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_splitclaw.png", dpi=130)
    return dict(err=float(err), alpha=np.round(alpha, 3).tolist(), sols=sols)


if __name__ == "__main__":
    print("claw", fig_claw())
    print("adjacent", fig_adjacent())
    print("splitclaw", fig_splitclaw())
