"""Pilot fix: trim each query to the region matched by a backbone HMM before adding.

trim_queries(backbone_aln, queries, wd, log) builds one profile HMM on the whole
backbone (hmmbuild), runs nhmmer-free `hmmsearch --domtblout` of the queries
against it, keeps for each query the span from the smallest envelope start to
the largest envelope end over its domain hits (plus PAD letters each side; an
overhang shorter than MINCUT is not cut) and
writes trimmed queries. The cut-off flanks are saved to flanks.json so that
restore() can put them back as unaligned letters (each in its own column at the
left/right edge of the alignment, lower case) once the trimmed queries have
been added. Queries with no hit are left untouched.

restore(..., realign=True) ("tfa" variant) first looks for flanks that are
homologous to each other: an all-vs-all nhmmer search of the flanks on each
side (E < 1e-5, hit to a different flank). Flanks with such a hit are aligned
to each other with MAFFT (--auto) and placed as one block of columns; the
others stay unaligned. This targets queries that carry a second region (e.g.
another gene/domain) that is absent from the backbone, which every HMM-based
adder leaves as unaligned insertion.
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
MINCUT = 50  # only cut an overhang this long: HMM envelopes often stop a few letters short


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
            lo = lo if lo >= MINCUT else 0
            hi = hi if len(s) - hi >= MINCUT else len(s)
        else:
            lo, hi = 0, len(s)
        out[n] = s[lo:hi]
        flanks[n] = (s[:lo], s[hi:])
    tq = os.path.join(wd, "queries.trimmed.fasta")
    fasta.write(out, tq)
    json.dump(flanks, open(os.path.join(wd, "flanks.json"), "w"))
    return tq


def homologous_flanks(fl, wd, tag, minlen=30, evalue=1e-5):
    """Names of flanks with a significant nhmmer hit to another flank."""
    fl = {n: s for n, s in fl.items() if len(s) >= minlen}
    if len(fl) < 2:
        return set()
    fa, tbl = os.path.join(wd, tag + ".fa"), os.path.join(wd, tag + ".tbl")
    fasta.write(fl, fa)
    subprocess.run([ENV + "/nhmmer", "--qformat", "fasta", "--dna", "--cpu", "4", "-E", str(evalue), "--tblout", tbl, fa, fa],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    keep = set()
    for line in open(tbl):
        if not line.startswith("#"):
            t = line.split()
            if t[0] != t[2] and float(t[12]) < evalue:  # target, query, E-value
                keep.update((t[0], t[2]))
    return keep


def side_block(fl, names, wd, tag, realign):
    """Columns for one side: aligned block of homologous flanks + singletons."""
    keep = homologous_flanks(fl, wd, tag) if realign else set()
    rows = {n: "" for n in names}
    if len(keep) >= 2:
        fa, out = os.path.join(wd, tag + ".keep.fa"), os.path.join(wd, tag + ".keep.aln")
        fasta.write({n: fl[n] for n in keep}, fa)
        with open(out, "w") as f:
            subprocess.run(["mafft", "--auto", "--thread", "4", fa], stdout=f, stderr=subprocess.DEVNULL, check=True)
        al = fasta.read(out)
        w = len(next(iter(al.values())))
        for n in names:
            rows[n] = al[n].upper() if n in al else "-" * w
    else:
        keep = set()
    single = [n for n in names if n not in keep and fl.get(n)]
    tot = sum(len(fl[n]) for n in single)
    off = 0
    for n in names:
        s = fl.get(n, "") if n in single else ""
        rows[n] += "-" * off + s.lower() + "-" * (tot - off - len(s))
        off += len(s)
    return rows, len(keep)


def restore(aln_path, wd, out_path, realign=False):
    """Put cut flanks back at the edges (unaligned, or realigned if homologous)."""
    aln = fasta.read(aln_path)
    flanks = json.load(open(os.path.join(wd, "flanks.json")))
    names = list(aln)
    left, nl = side_block({n: v[0] for n, v in flanks.items()}, names, wd, "left", realign)
    right, nr = side_block({n: v[1] for n, v in flanks.items()}, names, wd, "right", realign)
    fasta.write({n: left[n] + aln[n] + right[n] for n in names}, out_path)
    json.dump({"realigned_left": nl, "realigned_right": nr}, open(os.path.join(wd, "realign.json"), "w"))
    return out_path
