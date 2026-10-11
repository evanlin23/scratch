# AI-assisted (Claude), exploration code for CS581 project
"""Merge lane: for every (dataset, backbone size, B) run vote.py variants on MAGUS's fixed subsets and score.

    python3 runall.py [--once]

Backbones: size 200 #1-20 = the bank draw (MAGUS's own, inputs/backbones), #21-40 and sizes 100/400 = bbgen.py
(/opt/work/bbsize/bb/<DS>/s<size>). A condition runs once all its backbone files exist. Restartable: rows
already in results.jsonl are skipped. Rows -> cs581/bbsize/results/results.jsonl.
"""
import json, os, shutil, subprocess, sys, time

W = "/opt/work/bbsize"
REPO = "/home/user/scratch"
VOTE = REPO + "/cs581/gcmvote/code/vote.py"
RES = REPO + "/cs581/bbsize/results/results.jsonl"
sys.path.insert(0, REPO + "/cs581/code")
from gcmx.bbtool_bench import acc_ref  # noqa: E402

DS = ["BBA0101_b20", "BBA0067_b20", "SIMHIGH_R1", "SIMMOD_R1", "1000M2_R0_B20", "16S.M_R0_B20"]
SHORT = {d: d.split("_b20")[0].split("_R0_B20")[0] for d in DS}
VARIANTS = ["magus", "frac0.2", "frac0.3", "frac0.4", "frac0.5", "hard-bb"]
# (size, B) in priority order
CONDS = [(200, 10), (100, 10), (400, 10), (200, 40), (400, 5), (200, 5), (200, 20)]


def bb_files(d, size, B):
    if size == 200:
        out = []
        for i in range(1, B + 1):
            p = (os.path.join(W, "bank", d, "inputs", "backbones", "backbone_%d_mafft.txt" % i) if i <= 20
                 else os.path.join(W, "bb", d, "s200", "backbone_%d_mafft.txt" % i))
            out.append(p)
        return out
    return [os.path.join(W, "bb", d, "s%d" % size, "backbone_%d_mafft.txt" % i) for i in range(1, B + 1)]


def done_rows():
    if not os.path.exists(RES):
        return set()
    return {(r["dataset"], r["size"], r["B"], r["variant"]) for r in map(json.loads, open(RES))}


def edge_stats(vd):
    import numpy as np
    E = np.load(os.path.join(vd, "edges.npz"))
    k, n, w = E["k"], E["n"], E["w"].astype(float)
    f = k / n
    return {"mean_kn": float(f.mean()), "wmean_kn": float((w * f).sum() / w.sum()),
            "frac_unanimous": float((k == n).mean()), "frac_k1": float((k == 1).mean()),
            "mean_n": float(n.mean())}


def run_cond(d, size, B, done):
    rep = os.path.join(W, "bank", d)
    bbd = os.path.join(W, "cond", d, "s%d_B%d" % (size, B))
    shutil.rmtree(bbd, ignore_errors=True)
    os.makedirs(bbd)
    for i, p in enumerate(bb_files(d, size, B), 1):
        os.symlink(p, os.path.join(bbd, "backbone_%d_mafft.txt" % i))
    for v in VARIANTS:
        if (SHORT[d], size, B, v) in done:
            continue
        tag = "bbs_s%d_B%d_%s" % (size, B, v)
        t0 = time.time()
        r = subprocess.run([sys.executable, VOTE, rep, v, "--bb", bbd, "--tag", tag, "--threads", "4"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("FAIL", d, size, B, v, r.stdout[-500:], r.stderr[-1500:], flush=True)
            continue
        vd = os.path.join(rep, "vote", tag)
        info = json.loads(r.stdout.strip().splitlines()[-1])
        model = json.load(open(os.path.join(vd, "model.json")))
        s = acc_ref(os.path.join(rep, "true.fasta"), info["out"])
        row = {"dataset": SHORT[d], "size": size, "B": B, "variant": v, **s, "merge_wall": info["merge_wall"],
               "edges": model["edges"], "kept_edges": model["kept_edges"],
               "kept_weight_frac": model["kept_weight_frac"], "cutoff_by_n": model["cutoff_by_n"],
               "fit_bb": {k: model["fit"]["bb"][k] for k in ("pi", "p1", "p0", "conc1", "conc0")},
               "load1": round(os.getloadavg()[0], 2)}
        if v == "magus":
            row.update(edge_stats(vd))
        else:
            os.remove(os.path.join(vd, "edges.npz"))
        with open(RES, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps({k: row[k] for k in ("dataset", "size", "B", "variant", "avgErr", "merge_wall")}),
              flush=True)
        done.add((SHORT[d], size, B, v))


def main():
    once = "--once" in sys.argv
    while True:
        done = done_rows()
        todo = [(d, s, b) for d in DS for (s, b) in CONDS
                if any((SHORT[d], s, b, v) not in done for v in VARIANTS)]
        if not todo:
            print("RUNALL_DONE", flush=True)
            return
        ready = [c for c in todo if all(os.path.exists(p) for p in bb_files(*c))]
        if not ready:
            if once:
                return
            time.sleep(60)
            continue
        run_cond(*ready[0], done)


if __name__ == "__main__":
    main()
