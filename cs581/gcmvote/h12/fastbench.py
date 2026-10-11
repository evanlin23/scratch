# AI-assisted (Claude), exploration code for CS581 project
"""gcmx.bbtool_bench with MAGUS's graph built by the vectorized builder (--gcmx-fastgraph true).

Same arguments as bbtool_bench; only the --gcmx-fastgraph flag passed to gcmx.run_magus changes.
Whether that graph is the same is checked afterwards: gg.py's `linsi` variant re-merges the same subsets
and backbones with MAGUS's original builder and must reproduce magus.fasta."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))
from gcmx import bbtool_bench  # noqa: E402

_timed = bbtool_bench.timed


def timed(cmd, log, env=None):
    cmd = ["true" if a == "false" and cmd[i - 1] == "--gcmx-fastgraph" else a for i, a in enumerate(cmd)]
    return _timed(cmd, log, env)


bbtool_bench.timed = timed
if __name__ == "__main__":
    bbtool_bench.main()
