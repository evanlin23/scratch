# AI-assisted (Claude), exploration code for CS581 project
"""Helper h13: one fresh MAGUS draw (paper flags) on RNASim 10k R0, a replicate dir (as
gcmgen/code/fresh.sh builds it), and a fastgraph-vs-MAGUS-graph equivalence control. Restartable.

    python3 magus_h13.py WORK TRUE_ALIGNMENT K

  WORK/d0/state.json   rows: magus (fastgraph, end to end), merge-mafft-slow (MAGUS's own graph
                       builder, merge only on the same subsets/backbones), graph_identical
  WORK/rep             replicate dir for gg.py / vote.py
Each step records wall, CPU, peak summed process-tree RSS and max single-process RSS (memrun.py).
"""
import filecmp
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "code"))
sys.path.insert(0, HERE)
from gcmx import fasta  # noqa: E402
from gcmx.bbtool_bench import MERGE_FLAGS, acc_ref, backbone_stats  # noqa: E402
from gcmx.e2e_bench import magus_flags  # noqa: E402
from memrun import run  # noqa: E402

CODE = os.path.join(REPO, "code")
T = "4"
# Restart of 2026-10-11 (orchestrator's option a): the paper's flags plus --recurse false, reusing the killed
# recursive run's top-level decomposition and finished backbones (H13_REUSE=1 keeps WORK/d0/magus).
EXTRA = os.environ.get("H13_EXTRA", "").split()
REUSE = os.environ.get("H13_REUSE") == "1"


def main(work, src, k):
    w = os.path.join(work, "d0")
    os.makedirs(w, exist_ok=True)
    state_path = os.path.join(w, "state.json")
    row = json.load(open(state_path)) if os.path.exists(state_path) else {}
    ref = fasta.upper(fasta.read(src))
    true, unaligned = os.path.join(w, "true.fasta"), os.path.join(w, "unaligned.fasta")
    seqs = fasta.ungap(ref)
    fasta.write(ref, true)
    fasta.write(seqs, unaligned)
    row.update({"dataset": "RNASim10k_R0", "draw": 0, "threads": 4, "nproc": os.cpu_count(), "nseq": len(seqs),
                "nref": len(ref), "src": src, "k": k})

    def log_row(key, data):
        row[key] = data
        json.dump(row, open(state_path + ".tmp", "w"))
        os.replace(state_path + ".tmp", state_path)
        print(json.dumps({key: data}), flush=True)

    py = [sys.executable, "-m", "gcmx.run_magus"]
    inputs = os.path.join(w, "inputs")
    if "magus" not in row:
        for d in (("inputs",) if REUSE else ("magus", "inputs")):
            shutil.rmtree(os.path.join(w, d), ignore_errors=True)
        out = os.path.join(w, "magus.fasta")
        if os.path.exists(out):
            os.remove(out)
        m = run(py + ["--gcmx-fastgraph", "true", "-np", T, "-d", os.path.join(w, "magus"), "-i", unaligned,
                      "-o", out] + magus_flags(k) + EXTRA, os.path.join(w, "magus.log"), cwd=CODE)
        shutil.copytree(os.path.join(w, "magus", "subalignments"), os.path.join(inputs, "subalignments"))
        os.makedirs(os.path.join(inputs, "backbones"))
        os.makedirs(os.path.join(w, "graph_fast"), exist_ok=True)
        for f in os.listdir(os.path.join(w, "magus", "graph")):
            p = os.path.join(w, "magus", "graph", f)
            if f.startswith("backbone_") and f.endswith(("_unalign.txt", "_mafft.txt")):
                shutil.copy(p, os.path.join(inputs, "backbones"))
            elif os.path.isfile(p) and not f.startswith("backbone_"):
                shutil.copy(p, os.path.join(w, "graph_fast"))
        shutil.rmtree(os.path.join(w, "magus"), ignore_errors=True)
        log_row("magus", {**m, "fastgraph": True, "extra_flags": EXTRA, "reused_partial_run": REUSE, **acc_ref(true, out),
                          **backbone_stats(ref, os.path.join(inputs, "backbones"))})

    rep = os.path.join(work, "rep")
    if not os.path.exists(os.path.join(rep, "sets", "s0")):
        r = run([sys.executable, os.path.join(REPO, "protcons", "code", "pc.py"), "rep", w, rep],
                os.path.join(w, "rep.log"), cwd=CODE)
        print(json.dumps({"rep": r}), flush=True)

    if "merge-mafft-slow" not in row and os.environ.get("H13_SLOW", "1") == "1":
        d = os.path.join(w, "merge-slow")
        shutil.rmtree(d, ignore_errors=True)
        bb = os.path.join(d, "bb")
        os.makedirs(bb)
        for f in os.listdir(os.path.join(inputs, "backbones")):
            if f.endswith("_mafft.txt"):
                shutil.copy(os.path.join(inputs, "backbones", f), bb)
        out = os.path.join(w, "merge-mafft-slow.fasta")
        m = run(py + ["--gcmx-fastgraph", "false", "-np", T, "-d", os.path.join(d, "magus"),
                      "-s", os.path.join(inputs, "subalignments"), "-b", bb, "-o", out] + MERGE_FLAGS,
                os.path.join(w, "merge-mafft-slow.log"), cwd=CODE)
        g = os.path.join(d, "magus", "graph")
        same = {}
        for f in sorted(os.listdir(os.path.join(w, "graph_fast"))):
            if os.path.exists(os.path.join(g, f)):
                same[f] = filecmp.cmp(os.path.join(g, f), os.path.join(w, "graph_fast", f), shallow=False)
        same["output_alignment"] = filecmp.cmp(out, os.path.join(w, "magus.fasta"), shallow=False)
        shutil.rmtree(d, ignore_errors=True)
        log_row("merge-mafft-slow", {**m, **acc_ref(true, out), "identical_to_fastgraph": same})


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), sys.argv[3])
