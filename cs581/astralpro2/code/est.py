"""(4) Estimated gene trees: GDL families with branch lengths -> AliSim (GTR+G4, L bp) ->
FastTree; then ASTRAL-Pro3 / ASTRAL-multi / DISCO+ASTRAL on true and estimated trees.
usage: python est.py SETTING REP NGENES OUT.jsonl     (restartable per (setting, rep, kind, method))"""
import glob, json, math, os, random, re, subprocess, sys
from phylo import parse_newick, rf_error
from gdlsim import yule_tree
from gdl_len import sim_family

B = "/opt/mm/root/envs/gdl/bin/"
setting, rep, NG, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
W = "/opt/runs/ap2/est/%s_%d" % (setting, rep)
os.makedirs(W, exist_ok=True)
rng = random.Random(77 + rep)
SCALE, L = 0.6, 400   # tree height in subs/site; alignment length

if setting.startswith("S25"):
    nwk = yule_tree(26, seed=500 + rep, height=1.0)
    st = parse_newick(nwk)
    lamH = 0.93 if setting == "S25" else 2.0
    sig = 0.0 if setting == "S25" else 1.0
    lam, mu = {}, {}
    for v in range(len(st.parent)):
        lam[v] = lamH * (math.exp(rng.gauss(0, sig)) if sig else 1)
        mu[v] = lamH * (math.exp(rng.gauss(0, sig)) if sig else 1)
elif setting == "Rown":   # family-R own-rooting failure: a=b=0.02, c=0.6, lamT=2 (theory.md sec 5)
    nwk = "(((A:1,B:1)x:0.5,C:1.5)y:1,D:2.5);"
    st = parse_newick(nwk)
    rates = {"A": (0, math.log(50)), "B": (0, math.log(50)), "C": (0, math.log(1 / 0.6)), "x": (0, 0), "y": (2, 0), "D": (0, 0)}
    lam = {v: rates.get(st.label[v], (0, 0))[0] for v in range(len(st.parent))}
    mu = {v: rates.get(st.label[v], (0, 0))[1] for v in range(len(st.parent))}
    lam[st.root] = mu[st.root] = 0
true_sp = parse_newick(nwk)
H = 2.5 if setting == "Rown" else 1.0

gt_path = os.path.join(W, "true.nwk")
if not os.path.exists(gt_path):
    fams = []
    while len(fams) < NG:
        r = sim_family(st, lam, mu, rng)
        if r is None:
            continue
        g, cnt = r
        if len(cnt) < 4 or sum(cnt.values()) < 4:
            continue
        fams.append(g)
    open(gt_path, "w").write("\n".join(fams) + "\n")
fams = [l.strip() for l in open(gt_path) if ";" in l]

est_path = os.path.join(W, "est.nwk")
if not os.path.exists(est_path):
    ests = []
    for i, g in enumerate(fams):
        sc = re.sub(r":([0-9.eE-]+)", lambda m: ":%.6f" % (float(m.group(1)) * SCALE / H + 1e-6), g)
        tf = os.path.join(W, "g%d.nwk" % i); open(tf, "w").write(sc + "\n")
        pre = os.path.join(W, "a%d" % i)
        subprocess.run([B + "iqtree3", "--alisim", pre, "-t", tf, "-m", "GTR{1,2,1,1,2,1}+G4{1}", "--length", str(L),
                        "-seed", str(i + 1000 * rep), "-af", "fasta", "-redo"], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        fa = pre + ".fa"
        res = subprocess.run([B + "FastTree", "-nt", "-gtr", "-quiet", "-nopr", fa], check=True, capture_output=True, text=True)
        ests.append(res.stdout.strip())
        for f in glob.glob(pre + "*") + [tf]:
            os.remove(f)
    open(est_path + ".tmp", "w").write("\n".join(ests) + "\n")
    os.rename(est_path + ".tmp", est_path)

done = set()
if os.path.exists(out):
    done = {(r["setting"], r["rep"], r["kind"], r["method"]) for r in map(json.loads, open(out))}
for kind, path in [("true", gt_path), ("est", est_path)]:
    trees = [re.sub(r":[0-9.eE-]+", "", l.strip()) for l in open(path) if ";" in l]
    labs = sorted(set(re.findall(r"[A-Za-z0-9]+_\d+", "".join(trees))))
    mp = os.path.join(W, "map.txt"); open(mp, "w").write("".join("%s %s\n" % (l, l.rsplit("_", 1)[0]) for l in labs))
    gi = os.path.join(W, "in_%s.nwk" % kind); open(gi, "w").write("\n".join(trees) + "\n")
    for m in ["astral_pro3", "astral_multi", "disco_astral"]:
        if (setting, rep, kind, m) in done:
            continue
        o = os.path.join(W, "sp_%s_%s.nwk" % (kind, m))
        q = dict(check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if m == "astral_pro3":
            subprocess.run([B + "astral-pro3", "-a", mp, "-t", "1", "-o", o, gi], **q)
        elif m == "astral_multi":
            subprocess.run([B + "astral4", "-a", mp, "-t", "1", "-o", o, gi], **q)
        else:
            dd = os.path.join(W, "disco_%s.nwk" % kind)
            subprocess.run(["python3", "/opt/src/DISCO/disco.py", "-i", gi, "-o", dd, "-d", "_"], **q)
            subprocess.run([B + "astral4", "-a", mp, "-t", "1", "-o", o, dd], **q)
        est = parse_newick(open(o).read().strip())
        fn = rf_error(est, true_sp)
        rec = {"setting": setting, "rep": rep, "kind": kind, "method": m, "NG": NG, "fn": fn[0], "rf": list(fn)}
        open(out, "a").write(json.dumps(rec) + "\n")
        print(json.dumps(rec), flush=True)
