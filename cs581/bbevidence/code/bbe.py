"""What makes good GCM evidence? Merge-only experiments on cached MAGUS inputs.

    python3 bbe.py prep REP_DIR TRUE_ALIGNMENT          # normalise a replicate (uppercase, backbone seq sets)
    python3 bbe.py run REP_DIR VARIANT [VARIANT ...]    # build evidence, GCM merge, score -> REP_DIR/results.jsonl
    python3 bbe.py diag REP_DIR VARIANT [...]           # alignment-graph diagnostics -> REP_DIR/diag.jsonl

REP_DIR holds inputs/subalignments (25 L-INS-i subset alignments) and inputs/backbones (MAGUS's
10 L-INS-i backbones, backbone_N_mafft.txt). Variants are named (see VARIANTS / parse_variant):

  linsi            control: MAGUS's own backbones (must reproduce MAGUS's cached result)
  T                10 backbones realigned by tool T on MAGUS's own sequence sets (T in TOOLS)
  T@sK             10 backbones by T on new random sequence sets (seed K, 8 per subset)
  T@nN             N backbones by T: MAGUS's 10 sets, then new sets (seeds 1, 2, ...)
  A+B              union of the evidence of variants A and B (e.g. linsi+clustalo)
  A~5              first 5 backbones of A only
  A|gapT           A with columns of gap fraction > T masked (made insertion columns)
  A|consT          A with columns whose pairs agree with < T of the other backbones masked
Aligned backbones are cached per (tool, sequence set) in REP_DIR/aligned/<tool>/, with their time.
"""

import collections
import concurrent.futures
import json
import os
import random
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, "..", "..", "code"))
sys.path.insert(0, CODE)
from gcmx import fasta, score  # noqa: E402

THREADS = 4
MAFFT = os.environ.get("BBE_MAFFT", "")


def mafft_bin():
    global MAFFT
    if not MAFFT:
        from magus.configuration import Configs
        MAFFT = Configs.mafftPath
    return MAFFT


LINSI = ["--localpair", "--maxiterate", "1000", "--ep", "0.123", "--quiet", "--thread", str(THREADS), "--anysymbol"]
# tool name -> argv (input file appended). MAFFT entries use MAGUS's bundled binary.
TOOLS = {
    "linsi": ["MAFFT"] + LINSI,
    "clustalo": ["clustalo", "--threads", "1", "--force", "--outfmt", "fa", "-i"],
    "clustalo-iter": ["clustalo", "--threads", "1", "--force", "--outfmt", "fa", "--iter", "2", "-i"],
    "famsa": ["/opt/mm/root/envs/bio/bin/famsa", "-t", "1"],  # famsa IN OUT
    "muscle-super5": ["/usr/bin/muscle", "-threads", "1", "-super5"],  # muscle -super5 IN -output OUT
    "muscle": ["/usr/bin/muscle", "-threads", "1", "-align"],
    "fftns2": ["MAFFT", "--retree", "2", "--maxiterate", "0", "--quiet", "--thread", str(THREADS), "--anysymbol"],
    "fftnsi": ["MAFFT", "--retree", "2", "--maxiterate", "1000", "--quiet", "--thread", str(THREADS), "--anysymbol"],
    "mafft-auto": ["MAFFT", "--auto", "--quiet", "--thread", str(THREADS), "--anysymbol"],
    "ginsi": ["MAFFT", "--globalpair", "--maxiterate", "1000", "--quiet", "--thread", str(THREADS), "--anysymbol"],
    "linsi-noep": ["MAFFT", "--localpair", "--maxiterate", "1000", "--quiet", "--thread", str(THREADS), "--anysymbol"],
    "linsi-op2": ["MAFFT"] + LINSI + ["--op", "2.0"],
    "linsi-op3": ["MAFFT"] + LINSI + ["--op", "3.0"],
    "linsi-op4": ["MAFFT"] + LINSI + ["--op", "4.0"],
    "linsi-ul4": ["MAFFT"] + LINSI + ["--unalignlevel", "0.4"],
    "linsi-ul8": ["MAFFT"] + LINSI + ["--unalignlevel", "0.8"],
    "linsi-ul4-lexp": ["MAFFT"] + LINSI + ["--unalignlevel", "0.4", "--leavegappyregion"],
    "fftns2-op3": ["MAFFT", "--retree", "2", "--maxiterate", "0", "--op", "3.0", "--quiet", "--thread", str(THREADS),
                   "--anysymbol"],
}
MERGE_FLAGS = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false",
               "-f", "4"]


# ---------------------------------------------------------------- replicate preparation

def prep(rep, true_alignment):
    """Uppercase everything (MAGUS's map compares letters), write the reference and the
    unaligned backbone sequence sets (MAGUS's own sets, recovered from its backbones)."""
    fasta.write(fasta.upper(fasta.read(true_alignment)), os.path.join(rep, "true.fasta"))
    for d in ("subalignments", "backbones"):
        p = os.path.join(rep, "inputs", d)
        for f in os.listdir(p):
            fasta.write(fasta.upper(fasta.read(os.path.join(p, f))), os.path.join(p, f))
    sets = os.path.join(rep, "sets", "s0")
    os.makedirs(sets, exist_ok=True)
    for f in sorted(os.listdir(os.path.join(rep, "inputs", "backbones"))):
        n = f.split("_")[1]
        fasta.write(fasta.ungap(fasta.read(os.path.join(rep, "inputs", "backbones", f))),
                    os.path.join(sets, "backbone_{}.fa".format(n)))
    lin = os.path.join(rep, "aligned", "linsi", "s0")
    os.makedirs(lin, exist_ok=True)
    for f in os.listdir(os.path.join(rep, "inputs", "backbones")):
        shutil.copy(os.path.join(rep, "inputs", "backbones", f),
                    os.path.join(lin, "backbone_{}.fa".format(f.split("_")[1])))


def subsets(rep):
    d = os.path.join(rep, "inputs", "subalignments")
    files = sorted(os.listdir(d), key=lambda f: int(f.split("_")[-1].split(".")[0]))
    return [(f, fasta.read(os.path.join(d, f))) for f in files]


def new_sets(rep, seed, n=10, per_subset=8):
    """n new random backbone sequence sets (8 per subset, as MAGUS draws them)."""
    d = os.path.join(rep, "sets", "s{}".format(seed))
    if os.path.isdir(d) and len(os.listdir(d)) == n:
        return d
    os.makedirs(d, exist_ok=True)
    unal = fasta.ungap(fasta.read(os.path.join(rep, "true.fasta")))
    rng = random.Random(1000 + seed)
    subs = [list(s) for _, s in subsets(rep)]
    for b in range(n):
        taxa = [t for s in subs for t in rng.sample(s, min(per_subset, len(s)))]
        fasta.write({t: unal[t] for t in taxa}, os.path.join(d, "backbone_{}.fa".format(b + 1)))
    return d


def align_set(rep, tool, seed):
    """Align the 10 sequence sets of seed with tool, 4 at a time (as MAGUS schedules its backbone
    tasks); cached. Returns (dir, wall seconds, summed per-backbone seconds)."""
    out = os.path.join(rep, "aligned", tool, "s{}".format(seed))
    meta = os.path.join(rep, "aligned", tool, "s{}.json".format(seed))
    if os.path.exists(meta):
        m = json.load(open(meta))
        return out, m["wall"], m["sum"]
    if tool == "linsi" and seed == 0:
        return out, None, None
    src = os.path.join(rep, "sets", "s{}".format(seed)) if seed else os.path.join(rep, "sets", "s0")
    if seed and not os.path.isdir(src):
        new_sets(rep, seed)
    os.makedirs(out, exist_ok=True)
    files = sorted(f for f in os.listdir(src) if f.endswith(".fa"))

    def one(f):
        argv = [mafft_bin() if a == "MAFFT" else a for a in TOOLS[tool]]
        dst = os.path.join(out, f)
        start = time.time()
        if tool == "famsa":
            subprocess.run(argv + [os.path.join(src, f), dst], check=True, capture_output=True)
        elif tool.startswith("muscle"):
            subprocess.run(argv + [os.path.join(src, f), "-output", dst], check=True, capture_output=True)
        else:
            with open(dst + ".tmp", "w") as o:
                subprocess.run(argv + [os.path.join(src, f)], stdout=o, stderr=subprocess.DEVNULL, check=True)
            os.replace(dst + ".tmp", dst)
        fasta.write(fasta.upper(fasta.read(dst)), dst)
        return time.time() - start

    start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as pool:
        per = list(pool.map(one, files))
    m = {"wall": round(time.time() - start, 1), "sum": round(sum(per), 1)}
    json.dump(m, open(meta, "w"))
    return out, m["wall"], m["sum"]


# ---------------------------------------------------------------- masking

def mask_columns(aln, cols):
    """Turn columns into insertion columns: letters lowercase, gaps '.' (MAGUS skips both and
    does not advance its column index on them)."""
    cols = set(cols)
    out = {}
    for t, s in aln.items():
        out[t] = "".join((c.lower() if c != "-" else ".") if i in cols else c for i, c in enumerate(s))
    return out


def gap_mask(aln, t):
    rows = list(aln.values())
    n, L = len(rows), len(rows[0])
    arr = np.frombuffer("".join(rows).encode(), dtype=np.uint8).reshape(n, L)
    gapfrac = (arr == ord("-")).mean(axis=0)
    return [i for i in range(L) if gapfrac[i] > t]


def residue_columns(aln):
    """taxon -> array: alignment column of each residue."""
    return {t: np.array([i for i, c in enumerate(s) if c not in "-."], dtype=np.int32) for t, s in aln.items()}


def consistency_mask(alns, t):
    """For each backbone, per column: fraction of its residue pairs that every other backbone
    containing both residues also aligns (pairwise-agreement consistency, as in T-Coffee/M-Coffee
    column scores). Columns below t are masked. Pairs are sampled (<=200 per column)."""
    rcs = [residue_columns(a) for a in alns]
    rng = random.Random(0)
    masks = []
    for bi, a in enumerate(alns):
        taxa = list(a)
        L = len(a[taxa[0]])
        col_res = [[] for _ in range(L)]
        for t in taxa:
            for k, c in enumerate(rcs[bi][t]):
                col_res[c].append((t, k))
        bad = []
        for c in range(L):
            res = col_res[c]
            if len(res) < 2:
                continue
            pairs = [(res[i], res[j]) for i in range(len(res)) for j in range(i + 1, len(res))]
            if len(pairs) > 200:
                pairs = rng.sample(pairs, 200)
            agree = tot = 0
            for (x, i), (y, j) in pairs:
                for bj, rc in enumerate(rcs):
                    if bj == bi or x not in rc or y not in rc:
                        continue
                    tot += 1
                    agree += rc[x][i] == rc[y][j]
            if tot and agree / tot < t:
                bad.append(c)
        masks.append(bad)
    return masks


def agreement_mask(a, o, t):
    """Columns of a whose residue pairs are aligned by o (same sequences) in less than fraction t."""
    rco = residue_columns(o)
    rca = residue_columns(a)
    L = len(next(iter(a.values())))
    col_res = [[] for _ in range(L)]
    for x in a:
        for k, c in enumerate(rca[x]):
            col_res[c].append(rco[x][k])
    bad = []
    for c, oc in enumerate(col_res):
        n = len(oc)
        if n < 2:
            continue
        cnt = collections.Counter(oc)
        agree = sum(v * (v - 1) // 2 for v in cnt.values())
        if agree / (n * (n - 1) // 2) < t:
            bad.append(c)
    return bad


# ---------------------------------------------------------------- variants

def parse_variant(rep, name):
    """Returns (list of (label, alignment dict)), backbone wall seconds, summed seconds."""
    if "|" in name:
        base, how = name.rsplit("|", 1)
        files, wall, tot = parse_variant(rep, base)
        alns = [a for _, a in files]
        if how.startswith("gap"):
            t = float(how[3:])
            masks = [gap_mask(a, t) for a in alns]
        elif how.startswith("cons"):
            masks = consistency_mask(alns, float(how[4:]))
        elif how.startswith("agree"):  # agree-TOOL-T: columns whose pairs TOOL's alignment of the same set confirms
            _, other, t = how.split("-")
            ofiles, _, _ = parse_variant(rep, base.replace(base.split("@")[0].split("~")[0].split("^")[0], other, 1))
            masks = [agreement_mask(a, o, float(t)) for a, (_, o) in zip(alns, ofiles)]
        else:
            raise SystemExit("unknown mask " + how)
        return [(lab, mask_columns(a, m)) for (lab, a), m in zip(files, masks)], wall, tot
    if "+" in name:
        out, wall, tot = [], 0, 0
        for part in name.split("+"):
            f, w, s = parse_variant(rep, part)
            out += [(part + "_" + lab, a) for lab, a in f]
            wall, tot = wall + (w or 0), tot + (s or 0)
        return out, wall, tot
    if "^" in name:
        base, k = name.split("^")
        f, w, s = parse_variant(rep, base)
        k = int(k)
        return f[-k:], (w or 0) * k / len(f), (s or 0) * k / len(f)
    if "~" in name:
        base, k = name.split("~")
        f, w, s = parse_variant(rep, base)
        k = int(k)
        return f[:k], (w or 0) * k / len(f), (s or 0) * k / len(f)
    tool, seeds = name, [0]
    if "@s" in name:
        tool, s = name.split("@s")
        seeds = [int(s)]
    elif "@n" in name:
        tool, n = name.split("@n")
        seeds = list(range(int(n) // 10))
    out, wall, tot = [], 0, 0
    for seed in seeds:
        d, w, s = align_set(rep, tool, seed)
        wall, tot = wall + (w or 0), tot + (s or 0)
        for f in sorted(os.listdir(d), key=lambda f: int(f.split("_")[1].split(".")[0])):
            if f.endswith(".fa"):
                out.append(("s{}_{}".format(seed, f[:-3]), fasta.read(os.path.join(d, f))))
    return out, wall, tot


def safe(name):
    return name.replace("|", "_m_").replace("+", "_p_").replace("~", "_k_").replace("^", "_l_")


def merge(rep, name, files):
    vd = os.path.join(rep, "variants", safe(name))
    shutil.rmtree(vd, ignore_errors=True)
    bb = os.path.join(vd, "bb")
    os.makedirs(bb)
    for lab, a in files:
        fasta.write(a, os.path.join(bb, lab + ".txt"))
    out = os.path.join(vd, "out.fasta")
    start = time.time()
    with open(os.path.join(vd, "magus.log"), "w") as log:
        subprocess.run([sys.executable, "-m", "gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", str(THREADS),
                        "-d", os.path.join(vd, "work"), "-s", os.path.join(rep, "inputs", "subalignments"),
                        "-b", bb, "-o", out] + MERGE_FLAGS, cwd=CODE, stdout=log, stderr=subprocess.STDOUT,
                       check=True)
    wall = round(time.time() - start, 1)
    shutil.rmtree(os.path.join(vd, "work"), ignore_errors=True)
    return out, wall


def run(rep, names):
    res = os.path.join(rep, "results.jsonl")
    done = set()
    if os.path.exists(res):
        done = {json.loads(l)["variant"] for l in open(res)}
    for name in names:
        if name in done:
            continue
        files, bb_wall, bb_sum = parse_variant(rep, name)
        out, m_wall = merge(rep, name, files)
        s = score.fastsp(os.path.join(rep, "true.fasta"), out)
        row = {"rep": os.path.basename(rep.rstrip("/")), "variant": name, "nbb": len(files),
               "bb_wall": bb_wall, "bb_sum": bb_sum, "merge_wall": m_wall,
               **{k: s[k] for k in ("SPFN", "SPFP", "avgErr", "TC", "estHom", "refHom")}}
        row.update(cross_scores(rep, out))
        with open(res, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


# ---------------------------------------------------------------- pair-level diagnostics

class Rep:
    """Reference columns of every residue, subset membership, subset-column of every residue."""

    def __init__(self, rep):
        self.rep = rep
        ref = fasta.read(os.path.join(rep, "true.fasta"))
        self.refcol = residue_columns(ref)
        rows = list(ref.values())
        arr = np.frombuffer("".join(rows).encode(), dtype=np.uint8).reshape(len(rows), -1)
        self.ref_gapfrac = (arr == ord("-")).mean(axis=0)
        self.subset_of, self.subcol, self.subsets = {}, {}, []
        offset = 0
        for i, (f, a) in enumerate(subsets(rep)):
            self.subsets.append(a)
            rc = residue_columns(a)
            for t in a:
                self.subset_of[t] = i
                self.subcol[t] = rc[t] + offset  # global node id = offset + subset column
            offset += len(next(iter(a.values())))
        self.nnodes = offset
        # reference column composition of every node: node -> Counter(refcol)
        self.node_ref = collections.defaultdict(collections.Counter)
        self.node_size = np.zeros(self.nnodes, dtype=np.int64)
        for t, nodes in self.subcol.items():
            for node, rc in zip(nodes, self.refcol[t]):
                self.node_ref[node][rc] += 1
                self.node_size[node] += 1


def aln_pairs(R, aln, cross_only=True):
    """Counts of aligned residue pairs of an alignment (est), split true/false vs the reference,
    cross-subset only. Returns dict with est pairs, true pairs, and ref pairs among the same taxa."""
    rc = residue_columns(aln)
    taxa = list(aln)
    # group residues by est column: list of (subset, refcol)
    est_groups = collections.defaultdict(list)
    ref_groups = collections.defaultdict(list)
    for t in taxa:
        sub = R.subset_of[t]
        for ec, tc in zip(rc[t], R.refcol[t]):
            est_groups[ec].append((sub, tc))
            ref_groups[tc].append(sub)

    def cross(lst):
        n = len(lst)
        c = collections.Counter(lst)
        return n * (n - 1) // 2 - sum(k * (k - 1) // 2 for k in c.values())

    est = tp = 0
    for g in est_groups.values():
        est += cross([s for s, _ in g])
        byref = collections.defaultdict(list)
        for s, tc in g:
            byref[tc].append(s)
        tp += sum(cross(v) for v in byref.values())
    ref = sum(cross(v) for v in ref_groups.values())
    return {"est": est, "tp": tp, "ref": ref}


def cross_scores(rep, out):
    """SPFN/SPFP of the final alignment on cross-subset pairs only (the pairs GCM decides)."""
    R = REPS.setdefault(rep, Rep(rep))
    p = aln_pairs(R, fasta.upper(fasta.read(out)))
    return {"xSPFN": round(1 - p["tp"] / p["ref"], 5), "xSPFP": round(1 - p["tp"] / p["est"], 5)}


REPS = {}


def graph_edges(R, files, per_file=False):
    """MAGUS's alignment graph restricted to cross-subset edges, built from per-backbone residue
    pairs (vectorised). Returns (keys int64 = a * nnodes + b with a < b, weight, true weight,
    number of backbones contributing)."""
    N = R.nnodes
    acc_k, acc_w, acc_t, acc_n = [], [], [], []
    for _, a in files:
        rc = residue_columns(a)
        cols, nodes, refs, subs = [], [], [], []
        for t in a:
            cols.append(rc[t]); nodes.append(R.subcol[t]); refs.append(R.refcol[t])
            subs.append(np.full(len(rc[t]), R.subset_of[t], dtype=np.int32))
        cols, nodes, refs, subs = map(np.concatenate, (cols, nodes, refs, subs))
        order = np.argsort(cols, kind="stable")
        cols, nodes, refs, subs = cols[order], nodes[order], refs[order], subs[order]
        bounds = np.flatnonzero(np.diff(cols)) + 1
        ks, ts = [], []
        for lo, hi in zip(np.r_[0, bounds], np.r_[bounds, len(cols)]):
            if hi - lo < 2:
                continue
            i, j = np.triu_indices(hi - lo, 1)
            i, j = i + lo, j + lo
            keep = subs[i] != subs[j]
            i, j = i[keep], j[keep]
            na, nb = np.minimum(nodes[i], nodes[j]), np.maximum(nodes[i], nodes[j])
            ks.append(na.astype(np.int64) * N + nb)
            ts.append(refs[i] == refs[j])
        k = np.concatenate(ks); t = np.concatenate(ts)
        uk, inv = np.unique(k, return_inverse=True)
        acc_k.append(uk)
        acc_w.append(np.bincount(inv).astype(np.int64))
        acc_t.append(np.bincount(inv, weights=t).astype(np.int64))
        acc_n.append(np.ones(len(uk), dtype=np.int64))
    k = np.concatenate(acc_k)
    uk, inv = np.unique(k, return_inverse=True)
    w = np.bincount(inv, weights=np.concatenate(acc_w)).astype(np.int64)
    wt = np.bincount(inv, weights=np.concatenate(acc_t)).astype(np.int64)
    nbb = np.bincount(inv, weights=np.concatenate(acc_n)).astype(np.int64)
    return uk, w, wt, nbb


def final_node_cols(R, out):
    """node -> final alignment column."""
    rc = residue_columns(fasta.upper(fasta.read(out)))
    m = np.full(R.nnodes, -1, dtype=np.int64)
    for t, nodes in R.subcol.items():
        m[nodes] = rc[t]
    return m


def within_agreement(R, files):
    """Do the backbones repeat the subset alignments' within-subset errors? For backbone residue
    pairs (x, y) in the same subset: P(backbone aligns them | subset aligns them falsely) and
    P(... | subset aligns them correctly)."""
    agree_false = n_false = agree_true = n_true = 0
    for _, a in files:
        rc = residue_columns(a)
        groups = collections.defaultdict(list)  # (subset, subset column) -> backbone cols & ref cols
        for t in a:
            s = R.subset_of[t]
            for c, node, tc in zip(rc[t], R.subcol[t], R.refcol[t]):
                groups[node].append((c, tc))
        for g in groups.values():
            for i in range(len(g)):
                for j in range(i + 1, len(g)):
                    same_bb = g[i][0] == g[j][0]
                    if g[i][1] == g[j][1]:
                        n_true += 1
                        agree_true += same_bb
                    else:
                        n_false += 1
                        agree_false += same_bb
    return {"bb_keeps_subset_true": round(agree_true / max(n_true, 1), 4),
            "bb_keeps_subset_false": round(agree_false / max(n_false, 1), 4), "subset_false_pairs": n_false}


def diag(rep, names):
    R = REPS.setdefault(rep, Rep(rep))
    res = os.path.join(rep, "diag.jsonl")
    done = set()
    if os.path.exists(res):
        done = {json.loads(l)["variant"] for l in open(res)}
    for name in names:
        if name in done:
            continue
        files, _, _ = parse_variant(rep, name)
        row = {"rep": os.path.basename(rep.rstrip("/")), "variant": name}
        # backbone accuracy on cross-subset pairs (mean over backbones)
        ps = [aln_pairs(R, a) for _, a in files]
        row["bb_xprec"] = round(np.mean([p["tp"] / max(p["est"], 1) for p in ps]), 4)
        row["bb_xrec"] = round(np.mean([p["tp"] / max(p["ref"], 1) for p in ps]), 4)
        row["bb_xpairs"] = int(np.mean([p["est"] for p in ps]))
        # where are the backbones' false cross-subset pairs? gappiness of the backbone column and
        # of the residues' reference columns
        fp_gappy_bb = fp_tot = tp_gappy_bb = tp_tot = 0
        fp_gappy_ref = 0
        for _, a in files:
            rows = list(a.values())
            arr = np.frombuffer("".join(rows).encode(), dtype=np.uint8).reshape(len(rows), -1)
            gf = (arr == ord("-")).mean(axis=0)
            rc = residue_columns(a)
            cols = collections.defaultdict(list)
            for t in a:
                for c, tc in zip(rc[t], R.refcol[t]):
                    cols[c].append((R.subset_of[t], tc))
            for c, g in cols.items():
                by = collections.Counter(g)
                n = len(g)
                subc = collections.Counter(s for s, _ in g)
                cross = n * (n - 1) // 2 - sum(k * (k - 1) // 2 for k in subc.values())
                reft = collections.defaultdict(collections.Counter)
                for s, tc in g:
                    reft[tc][s] += 1
                tp = 0
                for tc, sc in reft.items():
                    m = sum(sc.values())
                    tp += m * (m - 1) // 2 - sum(k * (k - 1) // 2 for k in sc.values())
                fp = cross - tp
                fp_tot += fp
                tp_tot += tp
                if gf[c] > 0.5:
                    fp_gappy_bb += fp
                    tp_gappy_bb += tp
                # false pairs whose residues sit in gappy reference columns (> 50% gaps)
                if fp:
                    gref = sum(1 for _, tc in g if R.ref_gapfrac[tc] > 0.5) / n
                    fp_gappy_ref += fp * gref
        row["fp_in_gappy_bbcols"] = round(fp_gappy_bb / max(fp_tot, 1), 4)
        row["tp_in_gappy_bbcols"] = round(tp_gappy_bb / max(tp_tot, 1), 4)
        row["fp_in_gappy_refcols"] = round(fp_gappy_ref / max(fp_tot, 1), 4)
        # alignment graph
        keys, w, wt, nbb = graph_edges(R, files)
        N = R.nnodes
        ka, kb = keys // N, keys % N
        row["edges"] = len(keys)
        row["weight"] = int(w.sum())
        row["unit_prec"] = round(wt.sum() / w.sum(), 4)  # fraction of evidence units that are true pairs
        # edge-level truth: >= half of the residue pairs between the two full subset columns homologous
        purity = np.zeros(len(keys))
        for idx in range(len(keys)):
            ra, rb = R.node_ref[ka[idx]], R.node_ref[kb[idx]]
            if len(ra) > len(rb):
                ra, rb = rb, ra
            purity[idx] = sum(v * rb.get(c, 0) for c, v in ra.items()) / (R.node_size[ka[idx]] * R.node_size[kb[idx]])
        true_edge = purity >= 0.5
        row["edges_true"] = int(true_edge.sum())
        row["edges_false"] = int((~true_edge).sum())
        row["w_false_edges"] = int(w[~true_edge].sum())
        row["w_true_edges"] = int(w[true_edge].sum())
        row["mean_w_true_edge"] = round(w[true_edge].mean(), 2)
        row["mean_w_false_edge"] = round(w[~true_edge].mean(), 2) if (~true_edge).any() else 0
        for k in (1, 2, 5, 10):
            row["false_edges_w>={}".format(k)] = int(((~true_edge) & (w >= k)).sum())
            row["true_edges_w>={}".format(k)] = int((true_edge & (w >= k)).sum())
        row["false_edge_in_>=3bb"] = int(((~true_edge) & (nbb >= 3)).sum())
        row["true_edge_in_>=3bb"] = int((true_edge & (nbb >= 3)).sum())
        row["false_w_share_multi_bb"] = round(w[(~true_edge) & (nbb >= 2)].sum() / max(w[~true_edge].sum(), 1), 4)
        # false evidence units (wrong residue pairs) vs the edge they sit on
        row["false_units"] = int((w - wt).sum())
        row["false_units_on_false_edges"] = int((w - wt)[~true_edge].sum())
        # survival into the final alignment (needs the merged output)
        vd = os.path.join(rep, "variants", safe(name))
        out = os.path.join(vd, "out.fasta")
        if os.path.exists(out):
            fc = final_node_cols(R, out)
            merged = fc[ka] == fc[kb]
            row["false_edges_merged"] = int((merged & ~true_edge).sum())
            row["true_edges_merged"] = int((merged & true_edge).sum())
            row["false_w_merged"] = int(w[merged & ~true_edge].sum())
            row["true_w_merged"] = int(w[merged & true_edge].sum())
            row["false_w_survival"] = round(w[merged & ~true_edge].sum() / max(w[~true_edge].sum(), 1), 4)
            row["true_w_survival"] = round(w[merged & true_edge].sum() / max(w[true_edge].sum(), 1), 4)
            # merged node pairs that had NO direct evidence edge (transitive merges by MCL/trace)
        row.update(within_agreement(R, files))
        with open(res, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    cmd, rep = sys.argv[1], os.path.abspath(sys.argv[2])
    if cmd == "prep":
        prep(rep, sys.argv[3])
    elif cmd == "run":
        run(rep, sys.argv[3:])
    elif cmd == "diag":
        diag(rep, sys.argv[3:])
    elif cmd == "align":  # pre-align only: bbe.py align REP tool seed
        print(align_set(rep, sys.argv[3], int(sys.argv[4])))
