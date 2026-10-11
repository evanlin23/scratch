# AI-assisted (Claude), exploration code for CS581 project
"""Run vote.py variants on a replicate, score them (FastSP vs the reference), append to REP/results.jsonl.

    python3 run.py REP [--B 5|10|20] VARIANT [VARIANT ...]

--B 10 = MAGUS's own 10 L-INS-i backbones; 5 = backbones 1-5 of those; 20 = those 10 plus 10 new L-INS-i
backbones on new random 200-sequence sets (8 per subset, as MAGUS draws them; bbe.align_set seed 1).
Restartable: (variant, B) pairs already in results.jsonl are skipped.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, "..", "..", "code"))
sys.path.insert(0, CODE)
from gcmx.bbtool_bench import acc_ref  # noqa: E402


def bbdir(rep, B):
    src = os.path.join(rep, "inputs", "backbones")
    if B == 10:
        return src
    d = os.path.join(rep, "bb{}".format(B))
    if os.path.isdir(d) and len(os.listdir(d)) == B:
        return d
    orig = sorted(os.listdir(src), key=lambda f: int(f.split("_")[1]))
    tmp = d + ".tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    if B == 5:
        for f in orig[:5]:
            shutil.copy(os.path.join(src, f), tmp)
    elif B == 20:
        for f in orig:
            shutil.copy(os.path.join(src, f), tmp)
        sys.path.insert(0, os.path.join(HERE, "..", "..", "gcmgen", "code"))
        import gg  # noqa: F401  (patches bbe; pc.new_sets reads REP/unaligned.fasta when present)
        out, wall, tot = gg.bbe.align_set(rep, "linsi", 1)
        for f in sorted(os.listdir(out)):
            if f.endswith(".fa"):
                n = int(f.split("_")[1][:-3])
                shutil.copy(os.path.join(out, f), os.path.join(tmp, "backbone_{}_mafft.txt".format(10 + n)))
    else:
        raise SystemExit("B must be 5, 10 or 20")
    os.rename(tmp, d)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rep")
    ap.add_argument("variants", nargs="+")
    ap.add_argument("--B", type=int, default=10)
    a = ap.parse_args()
    rep = os.path.abspath(a.rep)
    res = os.path.join(rep, "results.jsonl")
    done = set()
    if os.path.exists(res):
        done = {(d["variant"], d.get("B")) for d in map(json.loads, open(res))}
    for v in a.variants:
        if (v, a.B) in done:
            continue
        tag = v if a.B == 10 else "{}@B{}".format(v, a.B)
        bb = bbdir(rep, a.B)
        r = subprocess.run([sys.executable, os.path.join(HERE, "vote.py"), rep, v, "--bb", bb, "--tag", tag],
                           capture_output=True, text=True)
        if r.returncode:
            print("FAILED", rep, tag, r.stderr[-2000:], flush=True)
            continue
        info = json.loads(r.stdout.strip().splitlines()[-1])
        vd = os.path.join(rep, "vote", tag)
        m = json.load(open(os.path.join(vd, "model.json")))
        s = acc_ref(os.path.join(rep, "true.fasta"), info["out"])
        row = {"rep": os.path.basename(rep), "variant": v, "B": a.B, "merge_wall": info["merge_wall"],
               "kept_edges": m["kept_edges"], "edges": m["edges"], "kept_weight_frac": m["kept_weight_frac"],
               "cutoff_by_n": m["cutoff_by_n"], **s}
        if "masked" in info:
            row.update(masked_cols=info["masked_cols"], cols=info["cols"])
        with open(res, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
