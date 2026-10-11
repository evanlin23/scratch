"""Score a test alignment against a BAliBASE 3 XML reference.

SP  = fraction of aligned residue pairs in reference core columns that the test reproduces (bali_score SP).
TC  = fraction of reference core columns reproduced in full (bali_score TC).
SPFN/SPFP over all reference homologies (FastSP definitions, full reference, case-insensitive).
Usage: python3 bbscore.py ref.xml test.fasta   -> prints SP TC SPFN SPFP
"""
import re, sys
from itertools import combinations


def read_xml(path):
    txt = open(path).read()
    names = re.findall(r"<seq-name>\s*(\S+)\s*</seq-name>", txt)
    seqs = [re.sub(r"\s", "", s) for s in re.findall(r"<seq-data>(.*?)</seq-data>", txt, re.S)]
    core = None
    for blk in re.findall(r"<column-score>(.*?)</column-score>", txt, re.S):
        if re.search(r"<colsco-name>\s*coreblock", blk):
            core = [int(x) for x in re.search(r"<colsco-data>(.*?)</colsco-data>", blk, re.S).group(1).split()]
    return dict(zip(names, seqs)), core


def read_fasta(path):
    d, n = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            n = line[1:].split()[0]; d[n] = []
        elif n:
            d[n].append(line)
    return {k: "".join(v) for k, v in d.items()}


def columns(aln, names):
    """List of columns; each column maps seq index -> residue index (ungapped)."""
    L = len(aln[names[0]])
    pos = [0] * len(names)
    cols = []
    for c in range(L):
        col = {}
        for i, n in enumerate(names):
            ch = aln[n][c]
            if ch not in "-.":
                col[i] = pos[i]; pos[i] += 1
        cols.append(col)
    return cols


def pairs(col):
    items = sorted(col.items())
    return {(a, ra, b, rb) for (a, ra), (b, rb) in combinations(items, 2)}


def score(ref_xml, test_fa):
    ref, core = read_xml(ref_xml)
    test = read_fasta(test_fa)
    names = list(ref)
    missing = [n for n in names if n not in test]
    if missing:
        raise SystemExit(f"missing in test: {missing[:3]}")
    for n in names:  # residues must agree
        a = re.sub(r"[-.]", "", ref[n]).upper(); b = re.sub(r"[-.]", "", test[n]).upper()
        if a != b:
            raise SystemExit(f"sequence mismatch for {n}")
    rcols = columns(ref, names)
    tcols = columns(test, names)
    # map (seq, residue) -> test column id
    where = {}
    for j, col in enumerate(tcols):
        for i, r in col.items():
            where[(i, r)] = j
    tpairs = set()
    for col in tcols:
        tpairs |= pairs(col)
    rpairs_all = set()
    core_pairs = core_hit = 0
    tc_n = tc_hit = 0
    for c, col in enumerate(rcols):
        p = pairs(col)
        rpairs_all |= p
        if core and core[c] == 1 and len(col) >= 2:
            core_pairs += len(p)
            hit = len(p & tpairs)
            core_hit += hit
            tc_n += 1
            # column reproduced if all residues of the ref column are in the same test column
            tc_hit += int(hit == len(p))
    shared = len(rpairs_all & tpairs)
    spfn = 1 - shared / len(rpairs_all) if rpairs_all else 0
    spfp = 1 - shared / len(tpairs) if tpairs else 0
    sp = core_hit / core_pairs if core_pairs else float("nan")
    tc = tc_hit / tc_n if tc_n else float("nan")
    return sp, tc, spfn, spfp


if __name__ == "__main__":
    print("%.4f %.4f %.4f %.4f" % score(sys.argv[1], sys.argv[2]))
