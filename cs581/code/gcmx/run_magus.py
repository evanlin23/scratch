"""Run MAGUS with an experimental GCM variant.

    python -m gcmx.run_magus [--gcmx-weight SCHEME] [--gcmx-alpha A] [--gcmx-order ORDER] <regular magus args>

New trace method: --graphtracemethod progdp (see progressive.py); pair it
with --graphclustermethod none since it does not use MCL clusters.

Unknown arguments are passed through to MAGUS unchanged. Use a fresh -d
working directory per run: MAGUS reuses graph/cluster files it finds there.
"""

import argparse
import sys

from . import progressive, weighting


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--gcmx-weight", default="count", choices=weighting.SCHEMES)
    parser.add_argument("--gcmx-alpha", type=float, default=1.0)
    parser.add_argument("--gcmx-order", default="upgma", choices=("upgma", "grow"))
    known, rest = parser.parse_known_args()

    weighting.install(known.gcmx_weight, known.gcmx_alpha)
    progressive.install(known.gcmx_order)

    from magus.main import main as magus_main
    sys.argv = [sys.argv[0]] + rest
    magus_main()


if __name__ == "__main__":
    main()
