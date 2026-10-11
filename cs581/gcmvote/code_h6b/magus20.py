# AI-assisted (Claude), exploration code for CS581 project
"""gcmx.bbtool_bench with the paper's MAGUS flags except -r 20 (20 backbones of 200 sequences).
    python3 magus20.py JOBFILE OUT.jsonl WORKDIR [bbtool_bench options]   (run from cs581/code)"""
import sys

from gcmx import bbtool_bench, e2e_bench

_flags = e2e_bench.magus_flags


def magus_flags(k):
    f = _flags(k)
    f[f.index("-r") + 1] = "20"
    return f


bbtool_bench.magus_flags = magus_flags
bbtool_bench.main()
