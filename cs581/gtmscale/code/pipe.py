"""DTM pipeline at scale, one dataset directory at a time (restartable: every step is
skipped when its output exists; timings/memory are kept in DIR/steps.json).

DIR must contain aln.fa (aligned) and true.tre.
Usage: python3 pipe.py DIR GUIDE MAXSUB [arms]
  GUIDE: ft (FastTree -gtr -gamma), ftfast (FastTree -fastest -noml, ME only),
         kmer (rapidnj on Mash-style k-mer distances of the unaligned sequences)
  arms (comma list, default all): full_ft, full_iq, gtm, blendft, polishft, cft, blendml, tm

Arms:
  gtm      GTM (convex mode) on guide + IQ-TREE subset trees
  blendft  GTM-Blend-FT: FastTree ML search started from the GTM tree, with every split
           of every subset tree as a partial (0/1/-) topological constraint. Subsets may
           interleave (blend) as long as each subset tree stays induced.
  polishft control: the same FastTree search from the GTM tree WITHOUT constraints
  cft      constrained FastTree from scratch (no GTM tree)
  blendml  GTM-Blend-ML from the prior pilot (python likelihood + RAxML-NG), time-capped
  full_ft, full_iq: full-data baselines (IQ-TREE capped by IQ_CAP seconds)
"""
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../gtm/code"))
from phylo import read_tree, fn_fp, bipartitions, is_induced  # noqa: E402
from sim import centroid_decomp  # noqa: E402
from blend import resolve  # noqa: E402

BIO = os.environ.get("BIOBIN", "/opt/mm/root/envs/bio/bin")
IQ = f"{BIO}/iqtree3"
FT = os.environ.get("FASTTREE", "FastTree")
GTM = os.environ.get("GTM_SRC", "/opt/gtmdata/GTM_src")
THREADS = int(os.environ.get("THREADS", "4"))
SUBFAST = os.environ.get("SUBFAST", "1") == "1"  # IQ-TREE -fast for subset trees (same as full IQ-TREE)
IQ_CAP = int(os.environ.get("IQ_CAP", "3600"))
FTOPT = ["-nt", "-gtr", "-gamma", "-quiet", "-nopr"]


def read_fasta(p):
    seqs, name = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


class Steps:
    def __init__(self, d):
        self.p = f"{d}/steps.json"
        self.s = json.load(open(self.p)) if os.path.exists(self.p) else {}

    def save(self):
        json.dump(self.s, open(self.p + ".tmp", "w"), indent=1)
        os.replace(self.p + ".tmp", self.p)


def timed(cmd, stdout=None, timeout=None):
    """Run cmd under GNU time; return (wall_s, maxrss_MB, returncode)."""
    tf = f"/tmp/gtmscale_time_{os.getpid()}_{time.time_ns()}"
    t0 = time.time()
    try:
        r = subprocess.run(["/usr/bin/time", "-f", "%M", "-o", tf, *cmd], stdout=stdout,
                           stderr=subprocess.DEVNULL, timeout=timeout)
        rc = r.returncode
    except subprocess.TimeoutExpired:
        rc = -9
    wall = time.time() - t0
    try:
        mem = int(open(tf).read().strip().splitlines()[-1]) / 1024
        os.remove(tf)
    except Exception:
        mem = float("nan")
    return wall, mem, rc


def kmer_nj(aln, out, k=8):
    """Mash-style distances on unaligned k-mer sets, then rapidnj."""
    import numpy as np
    from scipy import sparse
    seqs = read_fasta(aln)
    names = sorted(seqs)
    enc = {c: i for i, c in enumerate("ACGT")}
    rows, cols = [], []
    sizes = []
    for i, n in enumerate(names):
        s = seqs[n].upper().replace("-", "").replace("U", "T")
        codes = np.array([enc.get(c, -1) for c in s], dtype=np.int64)
        if len(codes) < k:
            ks = np.array([], dtype=np.int64)
        else:
            w = np.lib.stride_tricks.sliding_window_view(codes, k)
            w = w[(w >= 0).all(1)]
            ks = np.unique((w * (4 ** np.arange(k - 1, -1, -1))).sum(1))
        rows.append(np.full(len(ks), i))
        cols.append(ks)
        sizes.append(len(ks))
    X = sparse.csr_matrix((np.ones(sum(sizes), dtype=np.float32), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(len(names), 4 ** k))
    inter = (X @ X.T).toarray()
    sz = np.array(sizes, dtype=np.float64)
    union = sz[:, None] + sz[None, :] - inter
    jac = np.clip(inter / np.maximum(union, 1), 1e-6, 1)
    D = np.clip(-np.log(2 * jac / (1 + jac)) / k, 0, 5)
    np.fill_diagonal(D, 0)
    ph = out + ".phy"
    with open(ph, "w") as f:
        f.write("%d\n" % len(names))
        for i, n in enumerate(names):
            f.write(n + " " + " ".join("%.5f" % x for x in D[i]) + "\n")
    with open(out, "w") as f:
        subprocess.run([f"{BIO}/rapidnj", ph, "-i", "pd", "-n"], stdout=f, check=True, stderr=subprocess.DEVNULL)
    os.remove(ph)
    t = open(out).read().replace("'", "")
    open(out, "w").write(t)


def constraint_aln(subtrees, names, out):
    """One 0/1/- column per internal split of each subset tree."""
    idx = {n: i for i, n in enumerate(names)}
    cols = []
    for t in subtrees:
        leaves = [t.label[v] for v in t.leaves()]
        loc = {n: i for i, n in enumerate(leaves)}
        for b in bipartitions(t, loc):
            col = bytearray(b"-" * len(names))
            for n in leaves:
                col[idx[n]] = ord("1") if (b >> loc[n]) & 1 else ord("0")
            cols.append(col)
    with open(out, "w") as f:
        for i, n in enumerate(names):
            f.write(">%s\n%s\n" % (n, bytes(c[i] for c in cols).decode()))
    return len(cols)


def main():
    d, guide, maxsub = sys.argv[1], sys.argv[2], int(sys.argv[3])
    arms = (sys.argv[4] if len(sys.argv) > 4 else "full_ft,full_iq,gtm,blendft,polishft,cft").split(",")
    aln = f"{d}/aln.fa"
    if os.path.exists(f"{d}/SKIP"):  # budget: drop a replicate without editing running drivers
        return
    import fcntl
    lockf = open(f"{d}/.lock", "w")
    try:
        fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print("locked by another run:", d)
        return
    st = Steps(d)
    S = st.s
    T = read_tree(f"{d}/true.tre")

    def score(key, path):
        if os.path.exists(path) and "fn" not in S.get(key, {}):
            S[key]["fn"] = fn_fp(read_tree(path), T)[0]
            st.save()

    def run(key, cmd, out, stdout_to_out=False, timeout=None, extra=None):
        if key in S and os.path.exists(out):
            return
        if stdout_to_out:
            with open(out, "w") as f:
                w, m, rc = timed(cmd, stdout=f, timeout=timeout)
        else:
            w, m, rc = timed(cmd, timeout=timeout)
        S[key] = dict(wall=w, mem=m, rc=rc, **(extra or {}))
        if rc != 0 and rc != -9:
            raise RuntimeError(f"{key} failed rc={rc}: {cmd}")
        st.save()

    # full-data baselines (independent of guide)
    if "full_ft" in arms:
        run("full_ft", [FT, *FTOPT, aln], f"{d}/full_ft.tre", stdout_to_out=True)
        score("full_ft", f"{d}/full_ft.tre")
    if "full_iq" in arms and not os.path.exists(f"{d}/full_iq.treefile"):
        run("full_iq", [IQ, "-s", aln, "-m", "GTR+G", "-T", str(THREADS), "-fast", "-seed", "1", "--prefix",
                        f"{d}/full_iq", "-redo", "-quiet", "-mem", os.environ.get("IQ_MEM", "10G")], f"{d}/full_iq.treefile", timeout=IQ_CAP)
    if "full_iq" in arms:
        score("full_iq", f"{d}/full_iq.treefile")

    # guide
    g = f"{d}/guide_{guide}.tre"
    gk = f"guide_{guide}"
    if guide == "ft":
        if not os.path.exists(f"{d}/full_ft.tre"):
            run("full_ft", [FT, *FTOPT, aln], f"{d}/full_ft.tre", stdout_to_out=True)
        if not os.path.exists(g):
            os.symlink("full_ft.tre", g)
            S[gk] = dict(S["full_ft"])
            st.save()
    elif guide == "ftfast":
        run(gk, [FT, "-nt", "-fastest", "-noml", "-quiet", "-nopr", aln], g, stdout_to_out=True)
    elif guide.endswith("+it"):
        # iteration: the previous round's Blend-FT tree is the new guide (cost = previous round's total)
        base = guide[:-3]
        prev = f"{d}/{base}_m{maxsub}/blendft.tre"
        if not os.path.exists(g):
            os.symlink(os.path.relpath(prev, d), g)
            S[gk] = dict(wall=0.0, mem=float("nan"), rc=0)
            st.save()
    elif guide == "kmer":
        if gk not in S or not os.path.exists(g):
            t0 = time.time()
            kmer_nj(aln, g)
            S[gk] = dict(wall=time.time() - t0, mem=float("nan"), rc=0)
            st.save()
    score(gk, g)

    # decomposition + IQ-TREE subset trees (parallel, 1 thread each)
    w = f"{d}/{guide}_m{maxsub}"
    os.makedirs(w, exist_ok=True)
    GU = read_tree(g)
    sp = f"{w}/subsets.json"
    if not os.path.exists(sp):
        json.dump(centroid_decomp(GU, maxsub), open(sp, "w"))
    subsets = json.load(open(sp))
    seqs = read_fasta(aln)
    names = sorted(seqs)
    pre = f"{guide}_m{maxsub}_"

    def subtree(i):
        p = f"{w}/s{i}.tre"
        if os.path.exists(p):
            return None
        with open(f"{w}/s{i}.fa", "w") as f:
            for n in subsets[i]:
                f.write(">%s\n%s\n" % (n, seqs[n]))
        r = timed([IQ, "-s", f"{w}/s{i}.fa", "-m", "GTR+G", "-T", "1", "-seed", "1", "--prefix", f"{w}/s{i}",
                   "-quiet", "-redo", *(["-fast"] if SUBFAST else [])])
        os.replace(f"{w}/s{i}.treefile", p)
        for ext in ("iqtree", "log", "mldist", "ckp.gz", "bionj", "model.gz", "uniqueseq.phy"):
            if os.path.exists(f"{w}/s{i}.{ext}"):
                os.remove(f"{w}/s{i}.{ext}")
        return r
    if pre + "subtrees" not in S or not all(os.path.exists(f"{w}/s{i}.tre") for i in range(len(subsets))):
        t0 = time.time()
        with ThreadPoolExecutor(THREADS) as ex:
            rs = [r for r in ex.map(subtree, range(len(subsets))) if r]
        S[pre + "subtrees"] = dict(wall=time.time() - t0, mem=max([r[1] for r in rs] or [float("nan")]),
                                   k=len(subsets), sizes=[len(s) for s in subsets],
                                   cpu_serial=sum(r[0] for r in rs))
        st.save()
    subfiles = [f"{w}/s{i}.tre" for i in range(len(subsets))]
    subs = [read_tree(p) for p in subfiles]
    if "fn" not in S[pre + "subtrees"]:
        nfn = sum(fn_fp(s, T)[2] for s in subs)
        nint = sum(fn_fp(s, T)[4] for s in subs)
        S[pre + "subtrees"]["fn"] = nfn / nint
        st.save()

    # mergers
    gtm = f"{w}/gtm.tre"
    gtm_bin = f"{w}/gtm_bin.tre"  # GTM output with its few polytomies resolved (FastTree -intree needs degree 3)
    if "gtm" in arms or "blendft" in arms or "polishft" in arms or "blendml" in arms:
        run(pre + "gtm", [sys.executable, f"{GTM}/gtm.py", "-s", g, "-t", *subfiles, "-o", gtm], gtm)
        score(pre + "gtm", gtm)
        if not os.path.exists(gtm_bin):
            open(gtm_bin, "w").write(resolve(read_tree(gtm)).to_newick() + "\n")
    con = f"{w}/constraints.fa"
    if ("blendft" in arms or "cft" in arms or "blendfast" in arms) and not os.path.exists(con):
        t0 = time.time()
        nc = constraint_aln(subs, names, con)
        S[pre + "constraints"] = dict(wall=time.time() - t0, ncol=nc)
        st.save()
    if "blendft" in arms:
        run(pre + "blendft", [FT, *FTOPT, "-constraints", con, "-intree", gtm_bin, aln], f"{w}/blendft.tre",
            stdout_to_out=True)
        score(pre + "blendft", f"{w}/blendft.tre")
    if "blendfast" in arms:
        # cheaper blending: FastTree -fastest with 2 rounds of ML NNI, same constraints, from GTM
        run(pre + "blendfast", [FT, *FTOPT, "-fastest", "-mlnni", "2", "-constraints", con, "-intree", gtm_bin, aln],
            f"{w}/blendfast.tre", stdout_to_out=True)
        score(pre + "blendfast", f"{w}/blendfast.tre")
    if "polishft" in arms:
        run(pre + "polishft", [FT, *FTOPT, "-intree", gtm_bin, aln], f"{w}/polishft.tre", stdout_to_out=True)
        score(pre + "polishft", f"{w}/polishft.tre")
    if "cft" in arms:
        run(pre + "cft", [FT, *FTOPT, "-constraints", con, aln], f"{w}/cft.tre", stdout_to_out=True)
        score(pre + "cft", f"{w}/cft.tre")
    if "blendml" in arms:
        run(pre + "blendml", [sys.executable, f"{HERE}/run_blendml.py", aln, gtm_bin, f"{w}/blendml.tre", *subfiles],
            f"{w}/blendml.tre", timeout=int(os.environ.get("BLENDML_CAP", "3600")))
        score(pre + "blendml", f"{w}/blendml.tre")
    if "tm" in arms and pre + "treemerge" not in S:
        # TreeMerge (original code + PAUP*), topological distances on the guide tree as in Park et al. 2021
        import shutil
        TM_PY = os.environ.get("TM_PY", "/opt/mm/root/envs/tm27/bin/python")
        TM_DIR = os.environ.get("TM_DIR", "/opt/tools/treemerge/python")
        leafid = {GU.label[v]: v for v in GU.label}
        gnames = sorted(leafid)
        mat = f"{w}/guide_node.mat"
        with open(mat, "w") as f:
            f.write(f"{len(gnames)}\n")
            for a in gnames:
                dist = {leafid[a]: 0}
                stk = [leafid[a]]
                while stk:
                    v = stk.pop()
                    for u in GU.adj[v]:
                        if u not in dist:
                            dist[u] = dist[v] + 1
                            stk.append(u)
                f.write(a + " " + " ".join(str(dist[leafid[b]]) for b in gnames) + "\n")
        open(mat + "_taxlist", "w").write("\n".join(gnames) + "\n")
        work = f"{w}/tm_work"
        shutil.rmtree(work, ignore_errors=True)
        os.makedirs(work)
        env_ld = os.path.dirname(TM_PY) + "/../lib"
        os.environ["LD_LIBRARY_PATH"] = env_ld
        wall, mem, rc = timed([TM_PY, f"{TM_DIR}/treemerge.py", "-s", g, "-m", mat, "-x", mat + "_taxlist",
                               "-o", f"{w}/treemerge.tre", "-p", "/opt/tools/paup/paup", "-w", work, "-t", *subfiles])
        del os.environ["LD_LIBRARY_PATH"]
        shutil.rmtree(work, ignore_errors=True)
        os.remove(mat)
        S[pre + "treemerge"] = dict(wall=wall, mem=mem, rc=rc)
        st.save()
        score(pre + "treemerge", f"{w}/treemerge.tre")
    # constraint satisfaction of the blended outputs
    for a in ("blendft", "blendfast", "cft", "blendml", "polishft"):
        k = pre + a
        p = f"{w}/{a}.tre"
        if k in S and os.path.exists(p) and "induced" not in S[k]:
            E = read_tree(p)
            S[k]["induced"] = sum(is_induced(E, s) for s in subs)
            S[k]["k"] = len(subs)
            st.save()
    print(json.dumps({k: {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a != "sizes"}
                      for k, v in S.items()}), flush=True)


if __name__ == "__main__":
    main()
