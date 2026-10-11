"""Stream selected per-replicate files from the MAGUS paper's Results.zip.

    python fetch_published.py OUTDIR DATASET:REPS [DATASET:REPS ...] [--files f1,f2]

REPS is a range "0-9" or list "0,3,5". Files land in OUTDIR/DATASET/R<i>/.
Existing files are skipped (restartable). Uses the HTTP range reader from
cs581/code/gcmx/validate_published.py, so the 2.8 GB zip is never downloaded.
"""
import io
import os
import sys
import zipfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))
from gcmx.validate_published import HttpFile, URL  # noqa: E402

FILES = ["true_align.txt", "true_tree.tre", "gcm.txt", "gcm_slow.txt", "pasta_align.txt",
         "pasta_1_align.txt", "pasta_3_gcm_align.txt", "pasta_4_align.txt", "pasta.tre"]


def reps(spec):
    if "-" in spec:
        a, b = spec.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in spec.split(",")]


def main():
    args = sys.argv[1:]
    files = FILES
    if "--files" in args:
        i = args.index("--files")
        files = args[i + 1].split(",")
        args = args[:i] + args[i + 2:]
    out, specs = args[0], args[1:]
    z = zipfile.ZipFile(io.BufferedReader(HttpFile(URL), buffer_size=1 << 22))
    names = set(z.namelist())
    for spec in specs:
        ds, r = spec.split(":")
        for rep in reps(r):
            d = os.path.join(out, ds, "R%d" % rep)
            os.makedirs(d, exist_ok=True)
            for f in files:
                src = "Outputs/%s/R%d/%s" % (ds, rep, f)
                dst = os.path.join(d, f)
                if os.path.exists(dst) or src not in names:
                    continue
                with open(dst + ".part", "wb") as o:
                    o.write(z.read(src))
                os.rename(dst + ".part", dst)
            print(ds, rep, sorted(os.listdir(d)), flush=True)


if __name__ == "__main__":
    main()
