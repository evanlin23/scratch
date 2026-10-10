"""One-factor-at-a-time design around a ROSE-1000M2-like baseline; restartable.

    python driver.py SIMDIR [--reps 5] [--n 500] [--methods ...] [--cells base,clock0,...]

Each cell changes one factor from BASE. Replicate r uses seed 100+r in every cell,
so cells on the clock axis share the same birth-death topology and node heights.
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

BASE = {"shape": "bd", "sigma": 0.3, "height": 1.1, "indel": 0.004, "size": "geo:9"}
# indel volume (rate x mean length) held at 0.036 when the length distribution changes
CELLS = {
    "base": {},
    "clock0": {"sigma": 0.0},
    "clock0.7": {"sigma": 0.7},
    "clock1.2": {"sigma": 1.2},
    "clock2.0": {"sigma": 2.0},
    "pow1.7": {"size": "pow:1.7/200", "indel": 0.00536},
    "pow1.5": {"size": "pow:1.5/200", "indel": 0.00331},
    "balanced": {"shape": "balanced"},
    "caterpillar": {"shape": "caterpillar"},
    "caterpillarbdh": {"shape": "caterpillar-bdh"},  # first design, abandoned: spine branches ~0
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("simdir")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--methods", default="true,famsa,mafft-auto,magus,pasta")
    ap.add_argument("--cells", default=",".join(CELLS))
    ap.add_argument("--k", type=int, default=25)
    ap.add_argument("--bb", type=int, default=200)
    ap.add_argument("--rep-list", default=None)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    for r in ([int(x) for x in a.rep_list.split(",")] if a.rep_list else range(a.reps)):  # rep-major order: every cell gets rep 0 before any gets rep 1
        for c in a.cells.split(","):
            p = {**BASE, **CELLS[c]}
            d = os.path.join(a.simdir, "{}_n{}_r{}".format(c, a.n, r))
            try:  # claim the dataset so several workers share one queue (delete *.lock of unfinished ones to resume)
                os.makedirs(d + ".lock")
            except FileExistsError:
                continue
            if not os.path.exists(os.path.join(d, "true.fasta")):
                subprocess.run([sys.executable, os.path.join(HERE, "simulate.py"), d, "--n", str(a.n), "--seed", str(100 + r)]
                               + [x for k, v in p.items() for x in ("--" + k, str(v))], check=True)
            subprocess.run([sys.executable, os.path.join(HERE, "run_methods.py"), d, "--methods", a.methods, "--k", str(a.k), "--bb", str(a.bb), "--threads", str(a.threads)])


if __name__ == "__main__":
    main()
