"""Compare stopping rules for k on identical DecoDiPhy trajectories (output of ksel_run.py).

Rules (each maps a trajectory k=1..K -> chosen k):
  paper     : DecoDiPhy default. Walk k=1,2,..; stop at the first k with min(p) < 0.01 and return k-1;
              also stop (return k) if loss < 1e-10.
  bic(lam)  : argmin_k  n*log(loss_k/n) + lam * 3k * log(n)            (lam tuned on train trees)
  ratio(tau): smallest k with loss_{k+1}/loss_k > tau                    (tau tuned on train trees)
  adj       : like paper, but ALSO stop (return k-1) as soon as the k-solution places two queries
              on adjacent edges (sharing a node). Motivation: Theorem 1 - adjacent placements are never
              identifiable, and spurious extra placements are 'star' splits around nodes next to
              existing placements. Parameter-free.
  adj+bic   : min(adj, bic)
Metrics: exact-k accuracy, mean |k_hat - k|, Jaccard(placement edges at k_hat, true edges).
Paired test: Wilcoxon signed-rank on per-instance Jaccard differences vs paper.
"""
import json, glob, os, sys, collections, math
import numpy as np
from scipy.stats import wilcoxon, binomtest

D = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else None
TRAIN = {"bees", "birds-jarvis", "1kp", "beetles", "hemipteroid"}

MAXSEED = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
runs = []
for fn in glob.glob(os.path.join(D, "*.json")):
    r = json.load(open(fn))
    if r["seed"] <= MAXSEED:
        runs.append(r)
n_failed = len([f for f in glob.glob(os.path.join(D, "*.failed"))])


def adjacent(par, a, b):
    pa, pb = par.get(a), par.get(b)
    return pa == pb or pa == b or pb == a


def has_adj(r, rd):
    A = rd["anchors"]; par = r["parent"]
    return any(adjacent(par, A[i], A[j]) for i in range(len(A)) for j in range(i + 1, len(A)))


def rule_paper(r, thr=0.01):
    R = r["rounds"]
    for i, rd in enumerate(R):
        if min(rd["p"]) < thr:
            return max(rd["k"] - 1, 1)
        if rd["loss"] < 1e-10:
            return rd["k"]
    return R[-1]["k"]


def rule_adj(r, thr=0.01):
    R = r["rounds"]
    for rd in R:
        if min(rd["p"]) < thr or (rd["k"] > 1 and has_adj(r, rd)):
            return max(rd["k"] - 1, 1)
        if rd["loss"] < 1e-10:
            return rd["k"]
    return R[-1]["k"]


def rule_bic(r, lam):
    n = r["n"]; best = None
    for rd in r["rounds"]:
        s = n * math.log(max(rd["loss"], 1e-300) / n) + lam * 3 * rd["k"] * math.log(n)
        if best is None or s < best[0]:
            best = (s, rd["k"])
    return best[1]


def rule_ratio(r, tau):
    R = r["rounds"]
    for i in range(len(R) - 1):
        if R[i]["loss"] < 1e-10 or R[i + 1]["loss"] / max(R[i]["loss"], 1e-300) > tau:
            return R[i]["k"]
    return R[-1]["k"]


def step_features(r, prev, cur):
    par = r["parent"]
    A = cur["anchors"]
    anyadj = any(adjacent(par, A[i], A[j]) for i in range(len(A)) for j in range(i + 1, len(A)))
    return [math.log10(max(min(cur["p"]), 1e-6)),
            (prev["ybar"] - cur["ybar"]) / max(prev["ybar"], 1e-9),
            math.log10(max(cur["loss"], 1e-300) / max(prev["loss"], 1e-300)),
            float(anyadj)]


def fit_logit(train_runs):
    from sklearn.linear_model import LogisticRegression
    X, y = [], []
    for r in train_runs:
        R = r["rounds"]
        for i in range(1, len(R)):
            X.append(step_features(r, R[i - 1], R[i])); y.append(R[i]["k"] > r["k"])
    return LogisticRegression(C=1.0, max_iter=1000).fit(np.array(X), np.array(y))


def rule_logit(r, model, thr=0.5):
    """walk k upward; stop at the first step the classifier calls spurious, return k-1."""
    R = r["rounds"]
    for i in range(1, len(R)):
        if R[i - 1]["loss"] < 1e-10:
            return R[i - 1]["k"]
        if model.predict_proba(np.array([step_features(r, R[i - 1], R[i])]))[0, 1] > thr:
            return R[i - 1]["k"]
    return R[-1]["k"]


def jac(r, k):
    rd = next(x for x in r["rounds"] if x["k"] == k)
    A, T = set(rd["anchors"]), set(r["true_anchors"])
    return len(A & T) / len(A | T)


def evaluate(rule, subset):
    ks, js = [], []
    for r in subset:
        kh = rule(r)
        ks.append((kh, r["k"])); js.append(jac(r, kh))
    ks = np.array(ks)
    return dict(acc=float(np.mean(ks[:, 0] == ks[:, 1])), mae=float(np.mean(np.abs(ks[:, 0] - ks[:, 1]))),
                bias=float(np.mean(ks[:, 0] - ks[:, 1])), jac=float(np.mean(js)), js=np.array(js), n=len(subset))


train = [r for r in runs if r["tree"] in TRAIN]
test = [r for r in runs if r["tree"] not in TRAIN]
lines = [f"runs: {len(runs)} (train {len(train)}, test {len(test)}); seeds <= {MAXSEED}; "
         f"runs where the authors' solver crashed (excluded): {n_failed}"]
# tune on train (per noise level is NOT done: one global parameter, noise level unknown in practice)
lams = [0.1, 0.25, 0.5, 1, 2, 4, 8, 16]
taus = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
best_lam = max(lams, key=lambda l: evaluate(lambda r: rule_bic(r, l), train)["jac"])
best_tau = max(taus, key=lambda t: evaluate(lambda r: rule_ratio(r, t), train)["jac"])
thrs = [0.001, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1]
best_thr = max(thrs, key=lambda t: evaluate(lambda r: rule_paper(r, t), train)["jac"])
best_thr_adj = max(thrs, key=lambda t: evaluate(lambda r: rule_adj(r, t), train)["jac"])
lines.append(f"tuned on train trees ({sorted(TRAIN)}): bic lam={best_lam}, ratio tau={best_tau}, "
             f"paper min_p={best_thr}, adj min_p={best_thr_adj}")
logit = fit_logit(train)
lines.append("logistic stop-rule coefficients [log10 min p, rel. ybar drop, log10 loss ratio, adjacent]: "
             + str(np.round(logit.coef_[0], 2).tolist()) + f" intercept {logit.intercept_[0]:.2f}")
rules = {
    "paper(min_p=0.01)": rule_paper,
    f"bic(lam={best_lam})": lambda r: rule_bic(r, best_lam),
    f"ratio(tau={best_tau})": lambda r: rule_ratio(r, best_tau),
    f"paper(min_p tuned={best_thr})": lambda r: rule_paper(r, best_thr),
    "adj (paper+adjacency stop)": rule_adj,
    f"adj(min_p tuned={best_thr_adj})": lambda r: rule_adj(r, best_thr_adj),
    "learned (logistic, 4 features)": lambda r: rule_logit(r, logit),
    "oracle(true k)": lambda r: r["k"],
}
for setname, S in [("TEST trees", test), ("ALL trees", runs)]:
    for noise in ["noise0", "noise1", "noise2", "all"]:
        sub = [r for r in S if noise == "all" or r["noise"] == noise]
        if not sub:
            continue
        lines.append(f"\n### {setname}, {noise} (n={len(sub)})")
        lines.append("| rule | exact-k acc | mean abs(k_hat-k) | bias | mean Jaccard | dJac vs paper | Wilcoxon p | wins/losses |")
        lines.append("|---|---|---|---|---|---|---|---|")
        base = evaluate(rule_paper, sub)
        for name, f in rules.items():
            e = evaluate(f, sub)
            diff = e["js"] - base["js"]
            nz = diff[np.abs(diff) > 1e-12]
            p = wilcoxon(nz).pvalue if len(nz) >= 5 else float("nan")
            lines.append(f"| {name} | {e['acc']:.3f} | {e['mae']:.2f} | {e['bias']:+.2f} | {e['jac']:.3f} | "
                         f"{diff.mean():+.3f} | {p:.2g} | {(nz > 0).sum()}/{(nz < 0).sum()} |")
# per true k on test, all noise
lines.append("\n### TEST trees, noisy only (noise1+noise2): exact-k accuracy by true k")
sub = [r for r in test if r["noise"] != "noise0"]
lines.append("| true k | n | " + " | ".join(rules) + " |")
lines.append("|---|---|" + "---|" * len(rules))
for k in sorted({r["k"] for r in sub}):
    s2 = [r for r in sub if r["k"] == k]
    lines.append(f"| {k} | {len(s2)} | " + " | ".join(f"{evaluate(f, s2)['acc']:.2f}" for f in rules.values()) + " |")
# how often does the true solution itself have adjacent placements? (generator avoids it)
na = sum(any(adjacent(r["parent"], a, b) for i, a in enumerate(r["true_anchors"]) for b in r["true_anchors"][i + 1:]) for r in runs)
lines.append(f"\ntrue placements with an adjacent pair: {na}/{len(runs)} (the authors' generator excludes near placements)")
txt = "\n".join(lines)
print(txt)
if OUT:
    open(OUT, "w").write(txt + "\n")
