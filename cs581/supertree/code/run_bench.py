"""Restartable benchmark driver. One JSON line per (dataset, replicate, method) in results/<tag>.jsonl.

usage: python3 run_bench.py <tag> <workers> <methods,comma,sep> <case_list_file>
case_list_file: lines "<case_id>\t<source_trees>\t<true_tree>"
"""
import sys, os, json, time, traceback, warnings
from concurrent.futures import ProcessPoolExecutor
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_trees, fn_fp
import methods as M

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
RUNS = "/opt/runs/bench"


def method_fn(name, src, wd):
    """Returns (out_tree_path, info dict)."""
    out = os.path.join(wd, name + ".tre")
    info = {}
    if name == "astral3":
        w, m = M.astral3(src, out, os.path.join(wd, name), mem="8g")
    elif name == "astral4":
        w, m = M.astral4(src, out, os.path.join(wd, name))
    elif name == "tqmc":
        w, m = M.tqmc(src, out, os.path.join(wd, name))
    elif name == "scs":
        w, m = M.scs(src, out, os.path.join(wd, name))
    elif name == "mrlft":
        w, m = M.mrl_fasttree(src, out, os.path.join(wd, name))
    elif name.startswith("dc"):
        # dc:<guide>:<maxsize>:<subset method>   guide is another method's output in wd
        _, guide, ms, sm = name.split(":")
        gpath = os.path.join(wd, guide + ".tre")
        if not os.path.exists(gpath):
            raise RuntimeError("guide missing: " + gpath)
        info = M.dc_gtm(src, gpath, out, os.path.join(wd, name.replace(":", "_")), max_size=int(ms),
                        subset_method=sm, workers=1)
        w, m = info["wall"], max(info["peak_sub_mb"] or 0, info["peak_gtm_mb"] or 0)
    else:
        raise ValueError(name)
    info.update({"wall": w, "peak_mb": m})
    return out, info


def do_case(args):
    tag, case, src, true_path, meths = args
    wd = os.path.join(RUNS, tag, case)
    os.makedirs(wd, exist_ok=True)
    true = read_trees(true_path)[0]
    rows = []
    for name in meths:
        try:
            out, info = method_fn(name, src, wd)
            est = read_trees(out)[0]
            err = fn_fp(true, est)
            rows.append({"tag": tag, "case": case, "method": name, **err, **info})
        except Exception as e:
            rows.append({"tag": tag, "case": case, "method": name, "error": repr(e)[:300]})
            traceback.print_exc()
        with open(os.path.join(RES, f"{tag}.jsonl"), "a") as f:
            f.write(json.dumps(rows[-1]) + "\n")
        print(tag, case, name, rows[-1].get("rf"), rows[-1].get("wall"), flush=True)
    return rows


def main():
    tag, workers, meths, cases = sys.argv[1], int(sys.argv[2]), sys.argv[3].split(","), sys.argv[4]
    os.makedirs(RES, exist_ok=True)
    done = set()
    rp = os.path.join(RES, f"{tag}.jsonl")
    if os.path.exists(rp):
        for line in open(rp):
            r = json.loads(line)
            if "error" not in r:
                done.add((r["case"], r["method"]))
    jobs = []
    for line in open(cases):
        case, src, tp = line.rstrip("\n").split("\t")
        todo = [m for m in meths if (case, m) not in done]
        if todo:
            jobs.append((tag, case, src, tp, todo))
    print(len(jobs), "cases to run", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        list(ex.map(do_case, jobs))


if __name__ == "__main__":
    main()
