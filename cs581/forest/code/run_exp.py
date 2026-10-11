"""Experiment driver (restartable): python run_exp.py OUT.jsonl [conditions-name] [nproc]"""
import sys, os, json, tempfile, traceback, shutil, warnings
from multiprocessing import Pool
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline as P  # noqa: E402

KS = [100, 300, 1000, 3000, 10000, 100000]
CONDS = {
    # name: (n, regime, model, ks, reps, small_grid)
    "main": [
        (50, "U:0.005:0.05", "JC", [100, 200, 500, 1000, 2000, 5000], 20, False),  # Kim et al. short
        (100, "U:0.05:0.1", "JC", KS, 20, False),    # Kim et al. long
        (100, "U:0.1:0.4", "JC", KS, 20, False),     # deep (saturation)
        (100, "UH:2.0:1.0", "K2P", KS, 20, False),   # ultrametric h=2, lognormal(1) rates, K2P
    ],
    "large": [
        (500, "U:0.1:0.4", "JC", [300, 3000], 5, True),
        (500, "UH:2.0:1.0", "K2P", [300, 3000], 5, True),
    ],
}


def key(n, k, reg, model, rep):
    return f"{n}|{k}|{reg}|{model}|{rep}"


def work(args):
    n, k, reg, model, rep, small = args
    seed = (hash((n, k, reg, model)) % 100000) * 100 + rep if False else \
        int.from_bytes(f"{n}{k}{reg}{model}".encode(), "little") % 1000003 * 100 + rep
    w = tempfile.mkdtemp(prefix="fx_")
    try:
        r = P.run_replicate(n, k, reg, model, seed, w, small_grid=small)
        r["rep"] = rep
        r["key"] = key(n, k, reg, model, rep)
        return r
    except Exception as e:
        return {"key": key(n, k, reg, model, rep), "error": repr(e), "tb": traceback.format_exc()}
    finally:
        shutil.rmtree(w, ignore_errors=True)


def main():
    out = sys.argv[1]
    cname = sys.argv[2] if len(sys.argv) > 2 else "main"
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try:
                d = json.loads(line)
                if "error" not in d:
                    done.add(d["key"])
            except Exception:
                pass
    tasks = []
    # interleave conditions by replicate so partial runs cover every cell
    for rep in range(max(c[4] for c in CONDS[cname])):
        for (n, reg, model, ks, reps, small) in CONDS[cname]:
            if rep >= reps:
                continue
            for k in ks:
                if key(n, k, reg, model, rep) not in done:
                    tasks.append((n, k, reg, model, rep, small))
    print(f"{len(tasks)} tasks ({len(done)} done)", flush=True)
    with Pool(nproc) as pool, open(out, "a") as fh:
        for r in pool.imap_unordered(work, tasks):
            fh.write(json.dumps(r, default=float) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
