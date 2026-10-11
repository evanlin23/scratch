"""gcmx.run_magus plus (a) stage/task profiling and (b) support-pruned GCM graphs.

    MF_PROF=prof.json MF_ESK=K python mfrun.py <gcmx.run_magus args>

MF_PROF: write a JSON profile at exit: stage wall times (decomposition, graph build, MCL, trace, write)
  and every external command MAGUS runs (start, end, child CPU from wait4, kind = subset / backbone /
  decomp / mcl / other, by the output file name).
MF_ESK=K (K > 1, needs --gcmx-fastgraph true): keep an off-diagonal graph entry only if at least K
  backbone files put the two residues in one column (the `#esK` rule of claude/cs581-gcmgen).
"""
import atexit
import json
import os
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "code")))

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

from gcmx import fastgraph, run_magus  # noqa: E402

T0 = time.time()
PROF = {"stages": {}, "tasks": [], "graph_add": []}
LOCK = threading.Lock()


def _kind(cmd, dest):
    b = os.path.basename(dest or "")
    if b.startswith("backbone_"):
        return "backbone"
    if b.startswith("subalignment_"):
        return "subset"
    if "mcl" in cmd.split()[0]:
        return "mcl"
    if "decomposition" in (dest or "") or "decomposition" in cmd:
        return "decomp"
    return "other"


def install_profiling():
    from magus.tools import external_tools
    from magus.configuration import Configs
    import shutil

    def runCommand(**kwargs):
        command = kwargs["command"]
        dest = next(iter(kwargs.get("fileCopyMap", {}).values()), None)
        Configs.log("Running an external tool, command: {}".format(command))
        start = time.time()
        p = subprocess.Popen(command, shell=True, cwd=kwargs["workingDir"], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, universal_newlines=True)
        out = p.stdout.read()
        _, status, ru = os.wait4(p.pid, 0)
        p.returncode = os.waitstatus_to_exitcode(status)
        end = time.time()
        with LOCK:
            PROF["tasks"].append({"kind": _kind(command, dest), "out": dest, "start": round(start - T0, 2),
                                  "end": round(end - T0, 2), "cpu": round(ru.ru_utime + ru.ru_stime, 2)})
        if p.returncode != 0:
            Configs.error("Command encountered error: {}\n{}".format(command, out))
            raise subprocess.CalledProcessError(p.returncode, command)
        for srcPath, destPath in kwargs.get("fileCopyMap", {}).items():
            shutil.move(srcPath, destPath)
    external_tools.runCommand = runCommand

    def timed(mod, name, label):
        orig = getattr(mod, name)

        def wrapper(*a, **k):
            s = time.time()
            try:
                return orig(*a, **k)
            finally:
                PROF["stages"].setdefault(label, []).append([round(s - T0, 2), round(time.time() - T0, 2)])
        setattr(mod, name, wrapper)

    from magus.align import aligner
    from magus.align.merge import merger, alignment_graph
    from magus.align.merge.graph_build import graph_builder as gb
    timed(aligner, "decomposeSequences", "decomposition")
    timed(merger, "buildGraph", "buildGraph")       # includes waiting for subsets and backbones
    timed(merger, "clusterGraph", "cluster")
    timed(merger, "findTrace", "trace")
    timed(merger, "optimizeTrace", "optimize")
    timed(merger, "writeAlignment", "write")
    timed(gb, "buildMatrix", "buildMatrix")          # includes waiting for backbone tasks
    timed(gb, "addAlignmentFileToGraph", "graph_add_python")  # MAGUS's pure-Python per-backbone add
    timed(alignment_graph.AlignmentGraph, "writeGraphToFile", "writeGraph")
    timed(alignment_graph.AlignmentGraph, "initializeMatrix", "initMatrix")

    def dump():
        PROF["total"] = round(time.time() - T0, 2)
        with open(os.environ["MF_PROF"], "w") as f:
            json.dump(PROF, f)
    atexit.register(dump)


def pruned_buildMatrix(context):
    """fastgraph.buildMatrix (unweighted path) + support filter."""
    from magus.configuration import Configs
    from magus.tasks import task
    K = int(os.environ.get("MF_ESK", "1"))
    graph = context.graph
    n = graph.matrixSize
    files = [t.outputFile for t in task.asCompleted(context.backboneTasks)]
    for path in context.backbonePaths:
        if path not in files:
            files.append(path)
    total = sp.csr_matrix((n, n), dtype=np.int64)
    support = sp.csr_matrix((n, n), dtype=np.int64)
    for f in files:
        A = fastgraph._column_matrix(fastgraph._backbone_alignmap(context, f), n)
        prod = (A.T @ A).tocsr()
        total = total + prod
        support = support + (prod > 0).astype(np.int64)
    coo = total.tocoo()
    keep = np.ones(coo.nnz, dtype=bool)
    if K > 1:
        sup = support.tocsr()
        sv = np.asarray(sup[coo.row, coo.col]).ravel()
        keep &= (sv >= K) | (coo.row == coo.col)
    if Configs.graphBuildRestrict:
        sub = np.array([s for s, _ in graph.matSubPosMap])
        keep &= (sub[coo.row] != sub[coo.col]) | (coo.row == coo.col)
    removed = int((~keep).sum())
    total = sp.csr_matrix((coo.data[keep], (coo.row[keep], coo.col[keep])), shape=(n, n))
    total.indices = total.indices.astype(np.int32)
    graph.matrix = fastgraph.CSRGraph(total)
    Configs.log("[mfrun] graph {} entries from {} backbones, K={} removed {}".format(total.nnz, len(files), K, removed))


if __name__ == "__main__":
    if int(os.environ.get("MF_ESK", "1")) > 1:
        fastgraph.buildMatrix = pruned_buildMatrix  # run_magus installs fastgraph.buildMatrix by name
    if os.environ.get("MF_PROF"):
        import magus.main as mm
        orig_main = mm.main

        def main():
            install_profiling()  # after run_magus's own patches
            orig_main()
        mm.main = main
    run_magus.main()
