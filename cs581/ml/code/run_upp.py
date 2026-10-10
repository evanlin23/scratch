"""Estimated alignment for a fragmentary replicate with UPP (Nguyen et al. 2015; sepp 4.5.6, env ml2).

    python run_upp.py DATASET REP [THREADS]
reads $MLDATA/<ds>/R<rep>/unaligned.fasta (or derives it from true_align.fasta), writes upp.fasta
(UPP's masked alignment, insertion columns removed) and upp_time.json (wall/CPU seconds).
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402

UPP = "/opt/mm/root/envs/ml2/bin/run_upp.py"


def main(ds, rep, threads="1"):
    d = os.path.join(rt.MLDATA, ds, "R%s" % rep)
    un = os.path.join(d, "unaligned.fasta")
    if not os.path.exists(un):
        names, seqs = rt.read_fasta(os.path.join(d, "true_align.fasta"))
        rt.write_fasta(un, names, [s.replace("-", "").replace(".", "") for s in seqs])
    work = tempfile.mkdtemp(prefix="upp_%s_%s_" % (ds, rep))
    os.makedirs(os.path.join(work, "out"))
    env = dict(os.environ, PATH="/opt/mm/root/envs/ml2/bin:" + os.environ["PATH"], CONDA_PREFIX="/opt/mm/root/envs/ml2")
    t = time.time()
    with open(os.path.join(work, "run.log"), "w") as log:
        p = subprocess.Popen(["python3", UPP, "-s", un, "-o", "upp", "-d", os.path.join(work, "out"), "-x", str(threads), "-m", "dna",
                              "-seed", "1"], stdout=log, stderr=log, env=env, cwd=work)
        _, status, ru = os.wait4(p.pid, 0)
    if os.waitstatus_to_exitcode(status) != 0:
        raise SystemExit("UPP failed, see " + work)
    # UPP's own python children are counted in ru (they are waited for by the parent process)
    shutil.copy(os.path.join(work, "out", "upp_alignment_masked.fasta"), os.path.join(d, "upp.fasta"))
    json.dump({"seconds": round(time.time() - t, 1), "cpu_seconds": round(ru.ru_utime + ru.ru_stime, 1)},
              open(os.path.join(d, "upp_time.json"), "w"))
    print(d, open(os.path.join(d, "upp_time.json")).read())


if __name__ == "__main__":
    main(*sys.argv[1:])
