"""Fetch per-replicate inputs from the MAGUS paper's published Results.zip.

    python fetch.py DATASET REP [DATASET REP ...]

Streams members of Results.zip (Illinois Data Bank, doi:10.13012/B2IDB-2643961_V1)
with HTTP range requests and caches them under $MLDATA (default /opt/data/mlcache):
    $MLDATA/<DATASET>/R<rep>/{true_align,gcm,pasta_align}.fasta, pasta.tre, true_tree.tre,
    pasta_time.txt (PASTA wall time in seconds, from pastatre_timing.txt)
"""
import io
import os
import sys
import time
import urllib.request
import zipfile

URL = "https://databank.illinois.edu/datafiles/x13co/download"
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36"
MLDATA = os.environ.get("MLDATA", "/opt/data/mlcache")
MEMBERS = {
    "true_align.txt": "true_align.fasta",
    "gcm.txt": "gcm.fasta",
    "pasta_align.txt": "pasta_align.fasta",
    "pasta.tre": "pasta.tre",
    "true_tree.tre": "true_tree.tre",
    "pastatre_timing.txt": "pasta_time.txt",
    "gcmtxt_timing.txt": "gcm_time.txt",
}


class HttpFile(io.RawIOBase):
    """Seekable read-only file over HTTP range requests (enough for zipfile)."""

    def __init__(self, url):
        self.url, self.pos = url, 0
        self.size = int(self._get("bytes=0-0").headers["Content-Range"].split("/")[1])

    def _get(self, rng):
        for attempt in range(6):
            try:
                req = urllib.request.Request(self.url, headers={"User-Agent": UA, "Range": rng})
                return urllib.request.urlopen(req, timeout=120)
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


def fetch(z, dataset, rep):
    out = os.path.join(MLDATA, dataset, "R%s" % rep)
    os.makedirs(out, exist_ok=True)
    for member, local in MEMBERS.items():
        path = os.path.join(out, local)
        if os.path.exists(path):
            continue
        data = z.read("Outputs/%s/R%s/%s" % (dataset, rep, member))
        with open(path + ".tmp", "wb") as f:
            f.write(data)
        os.replace(path + ".tmp", path)
    return out


if __name__ == "__main__":
    z = zipfile.ZipFile(io.BufferedReader(HttpFile(URL), buffer_size=1 << 20))
    args = sys.argv[1:]
    for ds, rep in zip(args[::2], args[1::2]):
        print(fetch(z, ds, rep), flush=True)
