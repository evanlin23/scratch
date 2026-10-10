"""Pilot fix: trim each query to the region matched by a backbone HMM before adding.

trim_queries(backbone_aln, queries, wd, log) builds one profile HMM on the whole
backbone (hmmbuild), runs nhmmer-free `hmmsearch --domtblout` of the queries
against it, keeps for each query the span from the smallest envelope start to
the largest envelope end over its domain hits (plus PAD letters each side) and
writes trimmed queries. The cut-off flanks are saved to flanks.json so that
restore() can put them back as unaligned letters (each in its own column at the
left/right edge of the alignment, lower case) once the trimmed queries have
been added. Queries with no hit are left untouched.
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))
from gcmx import fasta  # noqa: E402

ENV = "/opt/mm/root/envs/add/bin"
PAD = 10


def trim_queries(bb, q, wd, log, mode="trim", pad=PAD):
    hmm = os.path.join(wd, "bb.hmm")
    with open(log, "a") as f:
        subprocess.run([ENV + "/hmmbuild", "--dna", "--cpu", "4", hmm, bb], stdout=f, stderr=f, check=True)
        dom = os.path.join(wd, "dom.txt")
        subprocess.run([ENV + "/hmmsearch", "--cpu", "4", "--domtblout", dom, "-E", "10", "--domE", "10",
                        "--max", hmm, q], stdout=subprocess.DEVNULL, stderr=f, check=True)
    span = {}
    for line in open(dom):
        if line.startswith("#"):
            continue
        t = line.split()
        name, a, b = t[0], int(t[19]), int(t[20])  # env from / env to (1-based)
        lo, hi = span.get(name, (a, b))
        span[name] = (min(lo, a), max(hi, b))
    seqs = fasta.read(q)
    out, flanks = {}, {}
    for n, s in seqs.items():
        if n in span:
            lo, hi = span[n]
            lo, hi = max(0, lo - 1 - pad), min(len(s), hi + pad)
        else:
            lo, hi = 0, len(s)
        out[n] = s[lo:hi]
        flanks[n] = (s[:lo], s[hi:])
    tq = os.path.join(wd, "queries.trimmed.fasta")
    fasta.write(out, tq)
    json.dump(flanks, open(os.path.join(wd, "flanks.json"), "w"))
    return tq


def restore(aln_path, wd, out_path):
    """Put cut flanks back: each flank letter gets its own column at the edges."""
    aln = fasta.read(aln_path)
    flanks = json.load(open(os.path.join(wd, "flanks.json")))
    L = sum(len(v[0]) for v in flanks.values())
    R = sum(len(v[1]) for v in flanks.values())
    res, lo, ro = {}, 0, 0
    for n, s in aln.items():
        left, right = flanks.get(n, ("", ""))
        res[n] = ("-" * lo + left.lower() + "-" * (L - lo - len(left)) + s
                  + "-" * ro + right.lower() + "-" * (R - ro - len(right)))
        lo += len(left)
        ro += len(right)
    fasta.write(res, out_path)
    return out_path
