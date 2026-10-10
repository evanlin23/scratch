"""Run MAGUS with an experimental GCM variant.

    python -m gcmx.run_magus [--gcmx-weight SCHEME] [--gcmx-alpha A] [--gcmx-order ORDER] <regular magus args>

New trace method: --graphtracemethod progdp (see progressive.py); pair it
with --graphclustermethod none since it does not use MCL clusters.

Unknown arguments are passed through to MAGUS unchanged. Use a fresh -d
working directory per run: MAGUS reuses graph/cluster files it finds there.
"""

import argparse
import os
import re
import subprocess
import sys

from . import fastgraph, progressive, weighting

# Backbone aligners: argv before the input file; MAGUS's own is
# mafft --localpair --maxiterate 1000 --ep 0.123 --quiet --thread N --anysymbol
BACKBONE_TOOLS = {
    "mafft": None,
    # one thread: MAGUS already runs one backbone task per core, and Clustal Omega's OpenMP threads
    # make it slower when several run at once (120 x 380 aa backbone: 13.9 s with --threads 4 vs 1.95 s)
    "clustalo": lambda t: ["clustalo", "--threads", "1", "--force", "--outfmt", "fa", "-i"],
    "mafft-auto": lambda t: ["mafft", "--auto", "--quiet", "--thread", str(t), "--anysymbol"],
    "linsi-noep": lambda t: ["mafft", "--localpair", "--maxiterate", "1000", "--quiet", "--thread", str(t), "--anysymbol"],
    "ginsi": lambda t: ["mafft", "--globalpair", "--maxiterate", "1000", "--quiet", "--thread", str(t), "--anysymbol"],
}


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--gcmx-weight", default="count", choices=weighting.SCHEMES)
    parser.add_argument("--gcmx-alpha", type=float, default=1.0)
    parser.add_argument("--gcmx-order", default="upgma", choices=("upgma", "grow"))
    parser.add_argument("--gcmx-refine", type=int, default=0, help="leave-one-out refinement rounds (progdp)")
    parser.add_argument("--gcmx-fastgraph", default="true", choices=("true", "false"),
                        help="vectorized graph construction (identical graph, much faster)")
    parser.add_argument("--gcmx-mclthreads", type=int, default=0,
                        help="run MCL with -te N threads (same clustering; MAGUS itself runs MCL single-threaded)")
    parser.add_argument("--gcmx-backbonetool", default="mafft", choices=sorted(BACKBONE_TOOLS),
                        help="aligner for the GCM backbones only (subsets stay MAFFT L-INS-i); mafft = MAGUS's own")
    known, rest = parser.parse_known_args()

    if known.gcmx_fastgraph == "true" and known.gcmx_weight == "count":
        fastgraph.install()
    weighting.install(known.gcmx_weight, known.gcmx_alpha)
    progressive.install(known.gcmx_order, known.gcmx_refine)
    if known.gcmx_mclthreads > 0:
        _thread_mcl(known.gcmx_mclthreads)
    if known.gcmx_backbonetool != "mafft":
        _backbone_tool(known.gcmx_backbonetool)

    from magus.main import main as magus_main
    sys.argv = [sys.argv[0]] + rest
    magus_main()


def _thread_mcl(threads):
    from magus.tools import external_tools
    original = external_tools.runMcl

    def runMcl(matrixPath, inflation, workingDir, outputPath):
        task = original(matrixPath, inflation, workingDir, outputPath)
        task.taskArgs["command"] += " -te {}".format(threads)
        task.json = task.toJson()  # MAGUS serializes tasks at construction; worker threads run the JSON copy
        return task
    external_tools.runMcl = runMcl


def _backbone_tool(name):
    """Align the GCM backbones (graph/backbone_N_mafft.txt) with another tool. Subset
    alignments also go through buildMafftAlignment, so only backbone outputs are rerouted."""
    from magus.tools import external_tools
    original = external_tools.buildMafftAlignment

    def buildMafftAlignment(inputPath, outputPath, subtablePath=None):
        task = original(inputPath, outputPath, subtablePath)
        if subtablePath is None and re.fullmatch(r"backbone_\d+_mafft\.txt", os.path.basename(outputPath)):
            from magus.configuration import Configs
            mafft = Configs.mafftPath if name != "clustalo" else None
            argv = BACKBONE_TOOLS[name](Configs.numCores)
            if mafft:
                argv[0] = mafft  # MAGUS's bundled MAFFT, same binary as the default backbones
            temp = os.path.join(os.path.dirname(outputPath), "temp_" + os.path.basename(outputPath))
            task.taskArgs["command"] = subprocess.list2cmdline(argv + [inputPath]) + " > " + temp
            task.json = task.toJson()  # see _thread_mcl
        return task
    external_tools.buildMafftAlignment = buildMafftAlignment


if __name__ == "__main__":
    main()
