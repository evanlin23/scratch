"""Validation layer 1: rescore the MAGUS paper's published alignments.

    python -m gcmx.validate_published OUT.jsonl [DATASET_PREFIX ...]

Streams Results.zip from the Illinois Data Bank (doi:10.13012/B2IDB-2643961_V1)
with HTTP range requests (no 2.8 GB download), and for every replicate scores
the published MAGUS and PASTA alignments against the reference alignment with
FastSP. Rows already in OUT.jsonl are skipped, so it can be resumed.

Files per replicate in Results.zip:
  gcm.txt               MAGUS (Fast)          gcm_slow.txt        MAGUS (Slow)
  pasta_align.txt       PASTA, 3 iterations   pasta_4_align.txt   PASTA, 4 iterations
  pasta_1_align.txt     PASTA, 1 iteration    pasta_3_gcm_align.txt  PASTA(3)+GCM
  true_align.txt / true_align_clean.txt   reference (clean = length-filtered)
"""

import io
import json
import os
import re
import sys
import tempfile
import time
import urllib.request
import zipfile

from . import fasta, score

URL = "https://databank.illinois.edu/datafiles/x13co/download"
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36"
METHODS = {
    "gcm.txt": "MAGUS(Fast)",
    "gcm_slow.txt": "MAGUS(Slow)",
    "pasta_align.txt": "PASTA(3)",
    "pasta_4_align.txt": "PASTA(4)",
    "pasta_1_align.txt": "PASTA(1)",
    "pasta_3_gcm_align.txt": "PASTA(3)+GCM",
}


class HttpFile(io.RawIOBase):
    def __init__(self, url):
        self.url, self.pos = url, 0
        self.size = int(self._get("bytes=0-0").headers["Content-Range"].split("/")[1])

    def _get(self, rng):
        for attempt in range(6):
            try:
                return urllib.request.urlopen(urllib.request.Request(self.url, headers={"User-Agent": UA, "Range": rng}),
                                              timeout=120)
            except OSError:
                if attempt == 5:
                    raise
                time.sleep(2 ** (attempt + 1))

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def readinto(self, b):
        if self.pos >= self.size:
            return 0
        data = self._get("bytes={}-{}".format(self.pos, min(self.size, self.pos + len(b)) - 1)).read()
        b[:len(data)] = data
        self.pos += len(data)
        return len(data)


def subsets_from_log(text):
    match = re.search(r"maximum number of subsets (\d+)", text)
    return int(match.group(1)) if match else None


def main():
    out_path, prefixes = sys.argv[1], sys.argv[2:]
    done = set()
    if os.path.exists(out_path):
        with open(out_path) as f:
            done = {(r["dataset"], r["rep"], r["file"]) for r in map(json.loads, f)}

    z = zipfile.ZipFile(io.BufferedReader(HttpFile(URL), buffer_size=1 << 22))
    names = set(z.namelist())
    reps = sorted({n.rsplit("/", 1)[0] for n in names if n.count("/") == 3 and n.endswith("/true_align.txt")})
    if prefixes:
        reps = [r for r in reps if any(r.startswith("Outputs/" + p) for p in prefixes)]

    for rep_path in reps:
        _, dataset, rep = rep_path.split("/")
        todo = [f for f in METHODS if rep_path + "/" + f in names and (dataset, rep, f) not in done]
        if not todo:
            continue
        with tempfile.TemporaryDirectory() as tmp:
            refs = {}
            for ref in ("true_align_clean.txt", "true_align.txt"):
                if rep_path + "/" + ref in names:
                    refs[ref] = z.extract(rep_path + "/" + ref, tmp)
            for f in todo:
                est = z.extract(rep_path + "/" + f, tmp)
                taxa = set(fasta.read(est))
                ref_name = next((r for r, p in refs.items() if set(fasta.read(p)) == taxa), None)
                row = {"dataset": dataset, "rep": rep, "file": f, "method": METHODS[f], "reference": ref_name}
                # gcm.txt -> gcmtxt_timing.txt ; pasta_4_align.txt -> pasta_4tre_timing.txt
                if f.startswith("gcm"):
                    timing = rep_path + "/" + f[:-4] + "txt_timing.txt"
                else:
                    timing = rep_path + "/" + f[:-len("_align.txt")] + "tre_timing.txt"
                if timing in names:
                    row["published_seconds"] = float(z.read(timing).decode().strip() or "nan")
                log = rep_path + "/log_" + f
                if f.startswith("gcm") and log in names:
                    row["subsets"] = subsets_from_log(z.read(log).decode(errors="replace"))
                if ref_name is None:
                    row["error"] = "no reference with matching taxa"
                else:
                    try:
                        row.update(score.fastsp(refs[ref_name], est))
                    except Exception as e:  # FastSP failure (e.g. memory) is recorded, not fatal
                        row["error"] = repr(e)[:300]
                os.remove(est)
                print(json.dumps(row), flush=True)
                with open(out_path, "a") as out:
                    out.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()
