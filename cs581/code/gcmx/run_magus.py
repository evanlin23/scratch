"""Run MAGUS with an experimental GCM variant.

    python -m gcmx.run_magus [--gcmx-weight SCHEME] [--gcmx-alpha A] [--gcmx-order ORDER] <regular magus args>

New trace method: --graphtracemethod progdp (see progressive.py); pair it
with --graphclustermethod none since it does not use MCL clusters.

Unknown arguments are passed through to MAGUS unchanged. Use a fresh -d
working directory per run: MAGUS reuses graph/cluster files it finds there.
"""

import argparse
import sys

from . import fastgraph, progressive, weighting


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--gcmx-weight", default="count", choices=weighting.SCHEMES)
    parser.add_argument("--gcmx-alpha", type=float, default=1.0)
    parser.add_argument("--gcmx-order", default="upgma", choices=("upgma", "grow"))
    parser.add_argument("--gcmx-refine", type=int, default=0, help="leave-one-out refinement rounds (progdp)")
    parser.add_argument("--gcmx-fastgraph", default="true", choices=("true", "false"),
                        help="vectorized graph construction (identical graph, much faster)")
    known, rest = parser.parse_known_args()

    if known.gcmx_fastgraph == "true" and known.gcmx_weight == "count":
        fastgraph.install()
    weighting.install(known.gcmx_weight, known.gcmx_alpha)
    progressive.install(known.gcmx_order, known.gcmx_refine)

    from magus.main import main as magus_main
    sys.argv = [sys.argv[0]] + rest
    magus_main()


if __name__ == "__main__":
    main()
