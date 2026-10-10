"""Run aligners + FastTree on one simulated dataset; restartable per method.

    python run_methods.py DATASET_DIR [--methods mafft-auto,famsa,magus,pasta,mafft-linsi,true]
        [--threads 4] [--k 25]

DATASET_DIR has true.fasta, unaligned.fasta, model.tre (from simulate.py).
Writes DATASET_DIR/results.json: per method wall/cpu seconds, FastSP SPFN/SPFP,
FastTree (-nt -gtr) tree FN/FP vs the model tree. "true" = FastTree on the true
alignment (alignment-error-free tree baseline). Alignments are kept gzipped.
Commands follow cs581/code/gcmx/e2e_bench.py (PASTA 1.8.3 defaults, 3 iterations;
MAGUS(Fast) with the paper's flags).
"""
import argparse
import gzip
import json
import os
import resource
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CODE = os.path.join(REPO, "cs581", "code")
sys.path.insert(0, CODE)
from gcmx import score  # noqa: E402  (FastSP wrapper; not modified)
from gcmx.e2e_bench import magus_flags, PASTA_PY, PASTA  # noqa: E402

BIO = "/opt/mm/root/envs/bio/bin"


def timed(cmd, log, cwd=None, stdout_file=None):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.time()
    with open(log, "w") as lf:
        out = open(stdout_file, "w") if stdout_file else lf
        proc = subprocess.run(cmd, cwd=cwd, stdout=out, stderr=lf)
        if stdout_file:
            out.close()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    if proc.returncode != 0:
        raise RuntimeError("failed ({}): {}".format(proc.returncode, " ".join(cmd)))
    return round(time.time() - start, 1), round(cpu, 1)


def flags_m(k, bb):
    f = magus_flags(k)
    f[f.index("-m") + 1] = str(bb)  # backbone size, scaled with n (paper: 200 at n=1000, K=25)
    return f


def align(method, d, unaligned, out, T, k, bb=200):
    w = os.path.join(d, "work_" + method)
    shutil.rmtree(w, ignore_errors=True)
    os.makedirs(w)
    log = os.path.join(d, method + ".log")
    if method == "mafft-auto":
        r = timed(["mafft", "--auto", "--thread", T, unaligned], log, stdout_file=out)
    elif method == "mafft-linsi":
        r = timed(["mafft", "--localpair", "--maxiterate", "1000", "--thread", T, unaligned], log, stdout_file=out)
    elif method == "famsa":
        r = timed([BIO + "/famsa", "-t", T, unaligned, out], log)
    elif method in ("magus", "magus-k25"):
        if method == "magus-k25":  # the paper's literal flags (K=25, 200-seq backbones), unscaled
            k, bb = 25, 200
        r = timed([sys.executable, "-m", "gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", T, "-d", w,
                   "-i", unaligned, "-o", out, "--graphbuildhmmextend", "false"] + flags_m(k, bb), log, cwd=CODE)
    elif method == "pasta":
        r = timed([PASTA_PY, PASTA, "-i", unaligned, "-o", w, "-d", "dna", "--num-cpus", T, "--iter-limit", "3",
                   "--temporaries", os.path.join(w, "tmp"), "-j", "pastajob"], log)
        alns = [os.path.join(w, p) for p in os.listdir(w)
                if p.startswith("pastajob.marker001.") and p.endswith(".aln") and "masked" not in p]
        if len(alns) != 1:
            raise RuntimeError("PASTA output not found")
        shutil.copy(alns[0], out)
    else:
        raise ValueError(method)
    shutil.rmtree(w, ignore_errors=True)
    return r


def bipartitions(path, taxa_ns):
    import dendropy
    t = dendropy.Tree.get(path=path, schema="newick", taxon_namespace=taxa_ns, preserve_underscores=True,
                          rooting="force-unrooted")
    t.encode_bipartitions()
    return {b.split_bitmask for b in t.bipartition_encoding if not b.is_trivial()}, t


def tree_error(est, true):
    import dendropy
    ns = dendropy.TaxonNamespace()
    tb, _ = bipartitions(true, ns)
    eb, _ = bipartitions(est, ns)
    return {"treeFN": round(len(tb - eb) / len(tb), 5), "treeFP": round(len(eb - tb) / max(1, len(eb)), 5)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("d")
    ap.add_argument("--methods", default="true,famsa,mafft-auto,magus,pasta")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--k", type=int, default=25)
    ap.add_argument("--bb", type=int, default=200)
    a = ap.parse_args()
    d = os.path.abspath(a.d)
    T = str(a.threads)
    res_path = os.path.join(d, "results.json")
    res = json.load(open(res_path)) if os.path.exists(res_path) else {}
    true, unaligned, model = (os.path.join(d, x) for x in ("true.fasta", "unaligned.fasta", "model.tre"))

    def save():
        with open(res_path + ".tmp", "w") as f:
            json.dump(res, f, indent=1)
        os.replace(res_path + ".tmp", res_path)

    for m in a.methods.split(","):
        if m in res and "treeFN" in res[m]:
            continue
        row = {}
        try:
            if m == "true":
                aln = true
            else:
                aln = os.path.join(d, m + ".fasta")
                wall, cpu = align(m, d, unaligned, aln, T, a.k, a.bb)
                s = score.fastsp(true, aln)
                row.update({"wall": wall, "cpu": cpu, **{x: s[x] for x in ("SPFN", "SPFP", "avgErr", "TC", "LenEst", "LenRef")}})
            tre = os.path.join(d, m + ".fasttree.tre")
            fwall, _ = timed(["FastTree", "-nt", "-gtr", "-quiet", aln], os.path.join(d, m + ".fasttree.log"), stdout_file=tre)
            row.update({"fasttree_wall": fwall, **tree_error(tre, model)})
            if m != "true":
                with open(aln, "rb") as f, gzip.open(aln + ".gz", "wb") as g:
                    shutil.copyfileobj(f, g)
                os.remove(aln)
        except Exception as e:  # record and move on; a rerun retries it
            row = {"error": str(e)[:300]}
        res[m] = row
        save()
        print(json.dumps({os.path.basename(d): {m: row}}), flush=True)


if __name__ == "__main__":
    main()
