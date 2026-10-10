"""Predicted-3Di evidence inside MAGUS on a BAliBASE RV100 set (merge-only, cached MAGUS inputs).

    cd cs581/code && python3 ../prost3di/code/rv100.py SET     (e.g. BBA0067)

Needs: /opt/p3d/rv100runs/SET/inputs/{subalignments,backbones} (cached MAGUS run, cs581-worker-1 branch),
/opt/p3d/rv100p/SET.{aa,3di}.fa and SET.db (predict3di.sh). Restartable; results in
/opt/p3d/rv100runs/SET/results.json. Conditions (subsets / backbones):
  merge-A  AA L-INS-i subsets / AA L-INS-i backbones (= MAGUS; reproduces the cached run)
  merge-B  AA subsets / 3Di L-INS-i backbones (same sequence sets)
  merge-C  AA subsets / union of A and B backbones (20)
  merge-D  3Di subsets / 3Di backbones
  merge-E  3Di subsets / union backbones
  full-linsi, full-linsi3di, full-fm   one alignment of the whole set (4 threads)
"""
import concurrent.futures, json, os, shutil, subprocess, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from bbscore import read_fasta
from gcmx import score

S = sys.argv[1]
W = f"/opt/p3d/rv100runs/{S}"
REF = f"/home/user/scratch/cs581/data/balibase_clean/RV100_{S}.fasta"
AA = f"/opt/p3d/rv100p/{S}.aa.fa"; TDI = f"/opt/p3d/rv100p/{S}.3di.fa"
MAT = "/opt/p3d/mat/mat3di.mafft"
MERGE = ["--graphclustermethod", "mcl", "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false", "-f", "4"]
RES = os.path.join(W, "results.json")
res = json.load(open(RES)) if os.path.exists(RES) else {}
aa, tdi = read_fasta(AA), read_fasta(TDI)


def save():
    json.dump(res, open(RES + ".tmp", "w"), indent=1); os.replace(RES + ".tmp", RES)


def timed(cmd, **kw):
    t0 = time.time(); r0 = os.times()
    subprocess.run(cmd, check=True, **kw)
    r1 = os.times()
    return time.time() - t0, (r1.children_user - r0.children_user) + (r1.children_system - r0.children_system)


def write(d, names, path):
    with open(path, "w") as f:
        for n in names:
            f.write(f">{n}\n{d[n]}\n")


def align3di(src, dst, threads=1, ep=True):
    """Re-align the sequences of alignment file src on predicted 3Di (MAFFT L-INS-i, 3Di matrix); write AA."""
    names = list(read_fasta(src))
    tmp = dst + ".3di.in"
    write(tdi, names, tmp)
    cmd = ["mafft", "--localpair", "--maxiterate", "1000", "--quiet", "--thread", str(threads), "--aamatrix", MAT]
    if ep:
        cmd += ["--ep", "0.123"]  # MAGUS's own backbone flags
    out = subprocess.run(cmd + [tmp], check=True, capture_output=True, text=True).stdout
    os.remove(tmp)
    al = {}; n = None
    for line in out.splitlines():
        if line.startswith(">"):
            n = line[1:].split()[0]; al[n] = []
        else:
            al[n].append(line.strip())
    with open(dst, "w") as f:
        for n in names:
            row = "".join(al[n]); it = iter(aa[n])
            f.write(f">{n}\n" + "".join("-" if c == "-" else next(it) for c in row) + "\n")


def realign_dir(src, dst, suffix):
    if os.path.exists(os.path.join(dst, "done.json")):
        return json.load(open(os.path.join(dst, "done.json")))
    shutil.rmtree(dst, ignore_errors=True); os.makedirs(dst)
    files = [f for f in os.listdir(src) if f.endswith(suffix)]
    t0 = time.time(); r0 = os.times()
    with concurrent.futures.ThreadPoolExecutor(4) as ex:
        list(ex.map(lambda f: align3di(os.path.join(src, f), os.path.join(dst, f)), files))
    r1 = os.times()
    info = {"wall": time.time() - t0, "cpu": (r1.children_user - r0.children_user) + (r1.children_system - r0.children_system)}
    json.dump(info, open(os.path.join(dst, "done.json"), "w"))
    return info


def acc(path):
    s = score.fastsp(REF, path)
    return {k: s[k] for k in ("SPFN", "SPFP", "avgErr", "TC")}


inp = os.path.join(W, "inputs")
res["bb3di_time"] = realign_dir(os.path.join(inp, "backbones"), os.path.join(W, "bb3di"), "_mafft.txt")
res["sub3di_time"] = realign_dir(os.path.join(inp, "subalignments"), os.path.join(W, "sub3di"), ".txt")
save()
conds = {"merge-A": ("subalignments", ["backbones"]), "merge-B": ("subalignments", ["bb3di"]),
         "merge-C": ("subalignments", ["backbones", "bb3di"]), "merge-D": ("sub3di", ["bb3di"]),
         "merge-E": ("sub3di", ["backbones", "bb3di"])}
for key, (sub, bbs) in conds.items():
    if key in res:
        continue
    d = os.path.join(W, key); shutil.rmtree(d, ignore_errors=True)
    bbdir = os.path.join(d, "bb"); os.makedirs(bbdir)
    for b in bbs:
        src = os.path.join(inp, b) if b == "backbones" else os.path.join(W, b)
        for f in os.listdir(src):
            if f.endswith("_mafft.txt"):
                shutil.copy(os.path.join(src, f), os.path.join(bbdir, b + "_" + f))
    subdir = os.path.join(inp, sub) if sub == "subalignments" else os.path.join(W, sub)
    subs = os.path.join(d, "subs"); os.makedirs(subs)
    for f in os.listdir(subdir):
        if f.endswith(".txt"):
            shutil.copy(os.path.join(subdir, f), subs)
    out = os.path.join(W, key + ".fasta")
    if os.path.exists(out):
        os.remove(out)  # MAGUS silently skips a run whose output exists
    wall, cpu = timed([sys.executable, "-m", "gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", "4", "-d",
                       os.path.join(d, "magus"), "-s", subs, "-b", bbdir, "-o", out] + MERGE,
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    res[key] = {"merge_wall": wall, "merge_cpu": cpu, **acc(out)}
    shutil.rmtree(d, ignore_errors=True); save(); print(key, res[key], flush=True)

for key in ["full-linsi3di", "full-fm", "full-linsi"]:
    if key in res:
        continue
    out = os.path.join(W, key + ".fasta")
    if key == "full-linsi":
        with open(out, "w") as f:
            wall, cpu = timed(["mafft", "--localpair", "--maxiterate", "1000", "--quiet", "--thread", "4", AA], stdout=f)
    elif key == "full-linsi3di":
        t0 = time.time(); r0 = os.times()
        align3di(AA, out, threads=4, ep=False)
        r1 = os.times(); wall = time.time() - t0
        cpu = (r1.children_user - r0.children_user) + (r1.children_system - r0.children_system)
    else:
        d = os.path.join(W, "fmwork"); shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
        src = f"/opt/p3d/rv100p/{S}.db"
        for f in os.listdir(src):
            if not f.startswith("db_ss_h"):
                shutil.copy(os.path.join(src, f), d)
        wall, cpu = timed(["/opt/mm/root/envs/fs/bin/foldmason", "structuremsa", os.path.join(d, "db"),
                           os.path.join(d, "out"), "--threads", "4"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.copy(os.path.join(d, "out_aa.fa"), out); shutil.rmtree(d)
    res[key] = {"wall": wall, "cpu": cpu, **acc(out)}
    save(); print(key, res[key], flush=True)
