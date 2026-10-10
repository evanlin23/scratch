"""Polish a given tree with RAxML-NG and score it (alternative 2: GTM + polish; large datasets).

    python polish.py ALN START_TREE TRUE_TREE OUT.jsonl TAG [--threads T] [--mode fast|full]

fast = RAxML-NG fast mode (simplified topology optimisation, KH stop rule) from START_TREE;
full = standard RAxML-NG search from START_TREE. Records wall and CPU seconds (all threads),
final log-likelihood, and FN/RF of the start and polished trees vs TRUE_TREE.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402
import treeerr  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    for x in ("aln", "start", "true", "out", "tag"):
        ap.add_argument(x)
    ap.add_argument("--threads", default="1")
    ap.add_argument("--mode", default="fast")
    a = ap.parse_args()
    clean = a.aln + ".clean.fasta" if not a.aln.endswith(".clean.fasta") else a.aln
    if not os.path.exists(clean):
        rt.clean_alignment(a.aln, clean + ".tmp")
        os.replace(clean + ".tmp", clean)
    work = tempfile.mkdtemp(prefix="polish_%s_" % a.tag)
    extra = ["--opt-topology", "simplified", "--stop-rule", "kh-mult"] if a.mode == "fast" else []
    cmd = [rt.RAXMLNG, "--search", "--msa", clean, "--model", "GTR+G", "--threads", a.threads, "--seed", "1",
           "--tree", a.start, "--prefix", os.path.join(work, "rx")] + extra
    t = time.time()
    with open(os.path.join(work, "rx.out"), "w") as log:
        p = subprocess.Popen(cmd, stdout=log, stderr=log)
        _, status, ru = os.wait4(p.pid, 0)
    if os.waitstatus_to_exitcode(status) != 0:
        raise SystemExit("RAxML-NG failed, see " + work)
    out_tree = os.path.join(os.path.dirname(a.start), "%s.polished.tre" % a.tag)
    os.replace(os.path.join(work, "rx.raxml.bestTree"), out_tree)
    lnl = None
    for line in open(os.path.join(work, "rx.raxml.log")):
        if line.startswith("Final LogLikelihood:"):
            lnl = float(line.split(":")[1])
    e0, e1 = treeerr.error(a.true, a.start), treeerr.error(a.true, out_tree)
    row = {"tag": a.tag, "mode": a.mode, "threads": a.threads, "seconds": round(time.time() - t, 1),
           "cpu_seconds": round(ru.ru_utime + ru.ru_stime, 1), "lnl_tool": lnl,
           "start_fn": e0["fn_rate"], "fn_rate": e1["fn_rate"], "rf_rate": e1["rf_rate"], "tree": out_tree}
    with open(a.out, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
