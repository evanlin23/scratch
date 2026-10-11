"""Paired statistics: Wilcoxon signed-rank (exact when possible), sign test, mean/median
paired difference and a percentile bootstrap 95% CI of the mean difference."""
import numpy as np
from scipy.stats import wilcoxon, binomtest


def paired(a, b, seed=0):
    """a, b: per-replicate values (same order). Difference d = b - a."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = b - a
    nz = d[d != 0]
    out = dict(n=len(d), mean_a=a.mean(), mean_b=b.mean(), mean_diff=d.mean(), median_diff=float(np.median(d)),
               n_b_better=int((d < 0).sum()), n_b_worse=int((d > 0).sum()), n_tie=int((d == 0).sum()))
    if len(nz) >= 1:
        try:
            out["wilcoxon_p"] = float(wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)
        except ValueError:
            out["wilcoxon_p"] = 1.0
        out["sign_p"] = float(binomtest(int((nz < 0).sum()), len(nz), 0.5).pvalue)
    else:
        out["wilcoxon_p"] = out["sign_p"] = 1.0
    rng = np.random.default_rng(seed)
    bs = [rng.choice(d, len(d)).mean() for _ in range(5000)]
    out["ci95"] = (float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)))
    return out


def fmt(o, scale=100):
    return ("n=%d  A %.2f%%  B %.2f%%  diff(B-A) mean %+.2f median %+.2f [95%% CI %+.2f, %+.2f]  "
            "B better/worse/tie %d/%d/%d  Wilcoxon p=%.3g  sign p=%.3g" %
            (o["n"], scale * o["mean_a"], scale * o["mean_b"], scale * o["mean_diff"], scale * o["median_diff"],
             scale * o["ci95"][0], scale * o["ci95"][1], o["n_b_better"], o["n_b_worse"], o["n_tie"],
             o["wilcoxon_p"], o["sign_p"]))
