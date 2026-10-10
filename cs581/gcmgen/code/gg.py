"""Generalising consensus-filtered GCM evidence across data types (merge-only, paired vs MAGUS).

    python3 gg.py rep  NAME REP_DIR TRUE [WORKER_REF]   # replicate from a cached MAGUS run (inputs.tar.xz)
    python3 gg.py run  REP_DIR VARIANT [...]            # merge-only variants -> REP_DIR/results.jsonl
    python3 gg.py agree REP_DIR TOOL [...]              # agreement of linsi with TOOL -> REP_DIR/agree.jsonl

Builds on cs581/bbevidence/code/bbe.py (variant grammar, masking, intersection, merge) and
cs581/protcons/code/pc.py (edge-support merge). New here:

  tools   einsi (E-INS-i), ginsi (G-INS-i), linsi-op3, linsi-sh (L-INS-i on the backbone sequences in a
          shuffled input order), linsi-rt (shuffled order AND a random guide tree via --treein), plus
          every bbe tool (fftns2, fftns2-op3, clustalo, ...)
  softW:A&B   soft intersection: pairs confirmed by B keep weight 1, unconfirmed pairs get weight W.
              MAGUS's edge weight is a residue-pair count summed over backbones, so this is A once plus
              (1/W - 1) copies of A&B; all weights are scaled by 1/W, which MCL ignores (checked with dupK).
  wsoftW:A&B  the same through per-file integer weights (run_wmerge.py): A at weight 1 plus A&B at weight
              1/W - 1; also wsoftW:A|cons0.7 (masked pairs of A down-weighted instead of deleted)
  dupK:A      K copies of A (scale control for soft weighting)
  A#esK       GCM graph keeps only cross-subset edges supported by >= K backbones (pc.merge_edgesup)
"""

import json
import os
import random
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "bbevidence", "code"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "protcons", "code"))
import bbe  # noqa: E402
import pc  # noqa: E402  (also patches bbe.new_sets)
import pairdiff  # noqa: E402
from gcmx import fasta  # noqa: E402
from gcmx.bbtool_bench import acc_ref  # noqa: E402

T = str(bbe.THREADS)
bbe.TOOLS.update({
    "einsi": ["MAFFT", "--genafpair", "--maxiterate", "1000", "--ep", "0", "--quiet", "--thread", T, "--anysymbol"],
    "linsi-sh": ["MAFFT"] + bbe.LINSI,
    "linsi-rt": ["MAFFT"] + bbe.LINSI,
})
CUSTOM = {"linsi-sh", "linsi-rt"}
_orig_align_set = bbe.align_set


def random_tree(n, rng):
    """MAFFT --treein format (newick2mafft.rb): one merge per line, '%5d %5d %10.5f %10.5f', 1-based
    cluster ids, smaller id first; the merged cluster keeps the smaller id."""
    cl, out = list(range(1, n + 1)), []
    while len(cl) > 1:
        a, b = sorted(rng.sample(cl, 2))
        cl.remove(b)
        out.append("%5d %5d %10.5f %10.5f" % (a, b, 0.1, 0.1))
    return "\n".join(out) + "\n"


def align_set(rep, tool, seed, only=None):
    """bbe.align_set for the shuffled-order / random-tree L-INS-i tools (same cache layout)."""
    if tool not in CUSTOM:
        return _orig_align_set(rep, tool, seed, only)
    out = os.path.join(rep, "aligned", tool, "s{}".format(seed))
    meta = os.path.join(rep, "aligned", tool, "s{}.json".format(seed))
    times_path = os.path.join(rep, "aligned", tool, "s{}.times.json".format(seed))
    if os.path.exists(meta) and not only:
        m = json.load(open(meta))
        return out, m["wall"], m["sum"]
    src = os.path.join(rep, "sets", "s{}".format(seed))
    files = sorted((f for f in os.listdir(src) if f.endswith(".fa")), key=lambda f: int(f.split("_")[1][:-3]))
    if only:
        files = [f for f in files if int(f.split("_")[1][:-3]) in only]
    times = json.load(open(times_path)) if os.path.exists(times_path) else {}
    os.makedirs(out, exist_ok=True)
    todo = [f for f in files if not os.path.exists(os.path.join(out, f))]

    def one(f):
        b = int(f.split("_")[1][:-3])
        rng = random.Random(7919 * b + seed)
        seqs = fasta.read(os.path.join(src, f))
        taxa = list(seqs)
        rng.shuffle(taxa)
        dst = os.path.join(out, f)
        tmp = "{}.{}".format(dst, os.getpid())
        fasta.write({t: seqs[t] for t in taxa}, tmp + ".in")
        argv = [bbe.mafft_bin() if a == "MAFFT" else a for a in bbe.TOOLS[tool]]
        if tool == "linsi-rt":
            open(tmp + ".tree", "w").write(random_tree(len(taxa), rng))
            argv += ["--treein", tmp + ".tree"]
        start = time.time()
        with open(tmp, "w") as o:
            subprocess.run(argv + [tmp + ".in"], stdout=o, stderr=subprocess.DEVNULL, check=True)
        fasta.write(fasta.upper(fasta.read(tmp)), tmp + "2")
        os.replace(tmp + "2", dst)
        for x in (tmp, tmp + ".in", tmp + ".tree"):
            if os.path.exists(x):
                os.remove(x)
        return f, round(time.time() - start, 1)

    import concurrent.futures
    start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=bbe.THREADS) as pool:
        for f, s in pool.map(one, todo):
            times[f] = s
            json.dump(times, open(times_path, "w"))
    wall = round(time.time() - start, 1) if len(todo) == len(files) else None
    total = round(sum(times.get(f, 0) for f in files), 1)
    if not only:
        json.dump({"wall": wall, "sum": total}, open(meta, "w"))
    return out, wall, total


bbe.align_set = align_set


def files_of(rep, name):
    """(list of (label, alignment)), backbone wall, backbone summed seconds."""
    if name.startswith("soft"):
        w, expr = name[4:].split(":", 1)
        m = round(1 / float(w)) - 1
        a, b = expr.split("&", 1)
        A, wa, sa = bbe.parse_variant(rep, a)
        I, wi, si = bbe.parse_variant(rep, expr)
        out = [("u_" + lab, x) for lab, x in A]
        for j in range(m):
            out += [("c{}_{}".format(j, lab), x) for lab, x in I]
        return out, wi, si
    if name.startswith("wsoft"):  # same as softW, via integer per-file weights (run_wmerge.py)
        w, expr = name[5:].split(":", 1)
        m = round(1 / float(w)) - 1
        A, wa, sa = bbe.parse_variant(rep, expr.replace("|", "&").split("&", 1)[0])
        I, wi, si = bbe.parse_variant(rep, expr)
        return [("u_" + lab, x) for lab, x in A] + [("c_" + lab, x) for lab, x in I], wi, si, \
            {"c_" + lab + ".txt": m for lab, _ in I}
    if name.startswith("dup"):
        k, expr = name[3:].split(":", 1)
        A, wa, sa = bbe.parse_variant(rep, expr)
        return [("d{}_{}".format(j, lab), x) for j in range(int(k)) for lab, x in A], wa, sa
    return bbe.parse_variant(rep, name)


def merge_weighted(rep, name, files, weights, k=1):
    """bbe.merge through run_wmerge.py (per-file integer weights)."""
    import shutil
    vd = os.path.join(rep, "variants", bbe.safe(name).replace(":", "_c_").replace("#", "_es_"))
    shutil.rmtree(vd, ignore_errors=True)
    bb = os.path.join(vd, "bb")
    os.makedirs(bb)
    for lab, a in files:
        fasta.write(a, os.path.join(bb, lab + ".txt"))
    json.dump(weights, open(os.path.join(vd, "weights.json"), "w"))
    out = os.path.join(vd, "out.fasta")
    start = time.time()
    with open(os.path.join(vd, "magus.log"), "w") as log:
        subprocess.run([sys.executable, os.path.join(HERE, "run_wmerge.py"), "--gcmx-fastgraph", "false",
                        "-np", str(bbe.THREADS), "-d", os.path.join(vd, "work"),
                        "-s", os.path.join(rep, "inputs", "subalignments"), "-b", bb, "-o", out] + bbe.MERGE_FLAGS,
                       cwd=bbe.CODE, stdout=log, stderr=subprocess.STDOUT, check=True,
                       env=dict(os.environ, GG_WEIGHTS=os.path.join(vd, "weights.json"), GG_ESK=str(k)))
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
        os.makedirs(os.path.join(rep, "variants"), exist_ok=True)
        lock = os.path.join(rep, "variants", bbe.safe(name).replace(":", "_c_").replace("#", "_es_") + ".lock")
        try:  # lanes may share a replicate: never run the same variant twice at once
            os.close(os.open(lock, os.O_CREAT | os.O_EXCL))
        except FileExistsError:
            continue
        start = time.time()
        base, k = (name.split("#es") + [None])[:2]
        got = files_of(rep, base)
        files, bb_wall, bb_sum = got[:3]
        prep_wall = round(time.time() - start, 1)
        if len(got) == 4:
            out, m_wall = merge_weighted(rep, name, files, got[3], int(k or 1))
        elif k:
            out, m_wall = pc.merge_edgesup(rep, name, files, int(k))
        else:
            out, m_wall = bbe.merge(rep, name.replace(":", "_c_"), files)
        s = acc_ref(os.path.join(rep, "true.fasta"), out)
        row = {"rep": os.path.basename(rep.rstrip("/")), "variant": name, "nbb": len(files), "bb_wall": bb_wall,
               "bb_sum": bb_sum, "prep_wall": prep_wall, "merge_wall": m_wall, **s}
        with open(res, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


def agree(rep, tools):
    """Reference-free agreement of MAGUS's L-INS-i backbones with each tool (cross-subset pairs):
    overlap = fraction of L-INS-i pairs the tool also aligns (mean over backbones, and pooled), the
    bbevidence support statistic of the unconfirmed pairs; plus, with the reference, the precision of
    the confirmed and the unconfirmed L-INS-i pairs and of the tool's own backbones (diagnostics only)."""
    path = os.path.join(rep, "agree.jsonl")
    done = {json.loads(l)["tool"] for l in open(path)} if os.path.exists(path) else set()
    R = bbe.Rep(rep)
    gid, off = {}, 0
    for t, rc in R.refcol.items():
        gid[t] = np.arange(off, off + len(rc), dtype=np.int64)
        off += len(rc)
    A, _, _ = bbe.parse_variant(rep, "linsi")
    PA = [pairdiff.pairs_of(R, a, gid) for _, a in A]
    for tool in tools:
        if tool in done:
            continue
        B, bw, bs = bbe.parse_variant(rep, tool)
        ov, acc = [], dict(a=0, a_tp=0, sh=0, sh_tp=0, ao=0, ao_tp=0, b=0, b_tp=0)
        for (ka, ta, _), (_, b) in zip(PA, B):
            kb, tb, _ = pairdiff.pairs_of(R, b, gid)
            inb = np.isin(ka, kb, assume_unique=True)
            ov.append(inb.mean())
            acc["a"] += len(ka); acc["a_tp"] += int(ta.sum())
            acc["sh"] += int(inb.sum()); acc["sh_tp"] += int(ta[inb].sum())
            acc["ao"] += int((~inb).sum()); acc["ao_tp"] += int(ta[~inb].sum())
            acc["b"] += len(kb); acc["b_tp"] += int(tb.sum())
        sup = pc.support(pc.FreeRep(rep), A, B)
        row = {"rep": os.path.basename(rep.rstrip("/")), "tool": tool, "overlap": round(float(np.mean(ov)), 4),
               "overlap_min": round(float(np.min(ov)), 4), "overlap_pooled": round(acc["sh"] / acc["a"], 4),
               "support_unconf": sup["support_a_only"], "support_conf": sup["support_shared"],
               "prec_linsi": round(acc["a_tp"] / acc["a"], 4), "prec_conf": round(acc["sh_tp"] / max(acc["sh"], 1), 4),
               "prec_unconf": round(acc["ao_tp"] / max(acc["ao"], 1), 4), "prec_tool": round(acc["b_tp"] / max(acc["b"], 1), 4),
               "pairs_tool_over_linsi": round(acc["b"] / acc["a"], 4), "bb_wall": bw, "bb_sum": bs,
               "bb_mean_s": round(bs / len(B), 1) if bs else None}
        with open(path, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


def make_rep(name, rep, true, ref=None):
    """Replicate from a cached MAGUS run on a claude/cs581-worker-* branch."""
    if os.path.exists(os.path.join(rep, "sets", "s0")):
        return
    git = lambda *a: subprocess.run(["git", "-C", HERE] + list(a), capture_output=True)
    refs = [ref] if ref else git("for-each-ref", "--format=%(refname:short)",
                                 "refs/remotes/origin/claude/cs581-worker-*").stdout.decode().split()
    for r in refs:
        blob = git("show", "{}:cs581/experiments/runs/{}/inputs.tar.xz".format(r, name))
        if blob.returncode == 0:
            break
    else:
        raise SystemExit("no cached inputs for " + name)
    os.makedirs(rep, exist_ok=True)
    subprocess.run(["tar", "xJ", "-C", rep], input=blob.stdout, check=True)
    prepj = git("show", "{}:cs581/experiments/runs/{}/prep.json".format(r, name)).stdout
    open(os.path.join(rep, "magus.json"), "wb").write(prepj)
    bbe.prep(rep, true)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "rep":
        make_rep(sys.argv[2], os.path.abspath(sys.argv[3]), sys.argv[4], *(sys.argv[5:6]))
    elif cmd == "run":
        run(os.path.abspath(sys.argv[2]), sys.argv[3:])
    elif cmd == "agree":
        agree(os.path.abspath(sys.argv[2]), sys.argv[3:])
