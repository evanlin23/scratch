r"""Search for counterexamples to (A) => (B) with a MILP over the competitor edge set S'.

For every unlabeled shape on n leaves and every matching S of size k WITHOUT a continuum
(L[V\V(S), I] has full column rank, i.e. condition (B)/(C) false), ask whether some S' with
|S'| <= k and some alpha give an exact alternative (condition (A)):
    g = -L alpha;  u notin V(S): g_u = 0 if u notin V(S'), g_u >= 1 if u in V(S')
                   u in V(S):    g_u <= -1 if u notin V(S'), free otherwise
    V(S') = endpoints of S' exactly, V(S') != V(S).
(strict inequalities are homogeneous, so ">= 1" is w.l.o.g. up to the big-M bound.)
A feasible MILP is a counterexample to the conjecture (A) <=> (B).
Also re-derives (B) <=> pattern classes for the table.

usage: python milp_enum.py NMIN NMAX KMIN KMAX OUT.md [lengths=unit|rand]
"""
import sys, os, itertools, collections, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "decodiphy", "code"))
from config_enum import shapes, laplacian

M = 1e3
ALLS = os.environ.get("ALLS") == "1"  # also non-matching S (node-measure identifiability)


def alt_milp(L, n, E, VS, k):
    N = L.shape[0]
    I = N - n
    m = len(E)
    G = -L[:, n:]
    # variable order: alpha (I), z_e (m), y_u (N)
    nv = I + m + N
    A, lo, hi = [], [], []

    def row():
        return np.zeros(nv)
    for u in range(N):
        r = row(); r[:I] = G[u]
        if u in VS:
            # g_u <= -1 + (M+1) y_u ; g_u >= -M
            r2 = r.copy(); r2[I + m + u] = -(M + 1); A.append(r2); lo.append(-np.inf); hi.append(-1)
            A.append(r.copy()); lo.append(-M); hi.append(np.inf)
        else:
            # g_u >= y_u ; g_u <= M y_u
            r2 = r.copy(); r2[I + m + u] = -1; A.append(r2); lo.append(0); hi.append(np.inf)
            r3 = r.copy(); r3[I + m + u] = -M; A.append(r3); lo.append(-np.inf); hi.append(0)
    adj = collections.defaultdict(list)
    for j, (a, b) in enumerate(E):
        adj[a].append(j); adj[b].append(j)
    for u in range(N):
        # y_u <= sum z_e ; y_u >= z_e
        r = row(); r[I + m + u] = 1
        for j in adj[u]:
            r[I + j] = -1
        A.append(r); lo.append(-np.inf); hi.append(0)
        for j in adj[u]:
            r = row(); r[I + m + u] = 1; r[I + j] = -1; A.append(r); lo.append(0); hi.append(np.inf)
    r = row(); r[I:I + m] = 1; A.append(r); lo.append(1); hi.append(k)
    # V(S') != V(S)
    r = row(); const = 0
    for u in range(N):
        if u in VS:
            r[I + m + u] = -1; const += 1
        else:
            r[I + m + u] = 1
    A.append(r); lo.append(1 - const); hi.append(np.inf)
    integrality = np.r_[np.zeros(I), np.ones(m + N)]
    bounds = Bounds(np.r_[-M * np.ones(I), np.zeros(m + N)], np.r_[M * np.ones(I), np.ones(m + N)])
    res = milp(np.zeros(nv), constraints=LinearConstraint(np.array(A), lo, hi), integrality=integrality,
               bounds=bounds, options={"time_limit": 60})
    if res.status == 0:
        z = res.x[I:I + m] > 0.5
        return [E[j] for j in range(m) if z[j]], res.x[:I]
    return None if res.status == 2 else "timeout"


def main():
    nmin, nmax, kmin, kmax, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    lt = sys.argv[6] if len(sys.argv) > 6 else "rand"
    rng = np.random.default_rng(7)
    lines = ["| n | k | shapes | matchings S | with continuum (B) | without continuum | (A) without (B): counterexamples | MILP timeouts | sec |",
             "|---|---|---|---|---|---|---|---|---|"]
    ce = []
    for n in range(nmin, nmax + 1):
        shp = shapes(n)
        for k in range(kmin, kmax + 1):
            t0 = time.time(); c = collections.Counter()
            for E in shp:
                N = 2 * n - 2
                lens = np.ones(len(E)) if lt == "unit" else rng.exponential(1, len(E)) + 0.05
                L = laplacian(N, E, lens)
                c["shapes"] += 1
                for S in itertools.combinations(range(len(E)), k):
                    VSl = [x for e in S for x in E[e]]
                    if len(set(VSl)) != len(VSl) and not ALLS:
                        continue
                    VS = set(VSl); c["S"] += 1
                    outside = [u for u in range(N) if u not in VS]
                    cont = np.linalg.matrix_rank(L[outside][:, n:]) < N - n
                    if cont:
                        c["cont"] += 1; continue
                    c["nocont"] += 1
                    r = alt_milp(L, n, E, VS, k)
                    if r == "timeout":
                        c["to"] += 1
                    elif r is not None:
                        c["ce"] += 1
                        if len(ce) < 20:
                            ce.append(f"n={n} k={k} E={E} S={[E[e] for e in S]} S'={r[0]} alpha={np.round(r[1], 3).tolist()}")
            lines.append(f"| {n} | {k} | {c['shapes']} | {c['S']} | {c['cont']} | {c['nocont']} | {c['ce']} | {c['to']} | {time.time() - t0:.0f} |")
            print(lines[-1], flush=True)
            open(out, "w").write("\n".join(lines) + "\n\ncounterexamples:\n" + "\n".join(ce) + "\n")


if __name__ == "__main__":
    main()
