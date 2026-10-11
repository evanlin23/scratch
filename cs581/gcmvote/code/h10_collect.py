"""AI-assisted (Claude), exploration code for CS581 project
Collect h10 alignment rows: gg baselines (W/reps/*/results.jsonl) + vote rows (W/vreps/*/results.jsonl), each tagged
with src, an order-independent md5 of its output alignment (aln_md5) and the other variants of the same replicate
whose output is the same alignment (identical_to).

    python3 h10_collect.py W OUT.jsonl
"""
import glob, hashlib, json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))
from gcmx import fasta  # noqa: E402

W, out = sys.argv[1:3]


def md5(path):
    if not os.path.exists(path):
        return None
    a = fasta.read(path)
    return hashlib.md5("".join(">{}\n{}\n".format(k, a[k].upper()) for k in sorted(a)).encode()).hexdigest()


rows = []
for f in sorted(glob.glob(os.path.join(W, "reps", "*", "results.jsonl"))):
    rep = os.path.dirname(f)
    for r in map(json.loads, open(f)):
        d = r["variant"].replace("#", "_es_")
        rows.append({"src": "gg", **r, "aln_md5": md5(os.path.join(rep, "variants", d, "out.fasta"))})
for f in sorted(glob.glob(os.path.join(W, "vreps", "*", "results.jsonl"))):
    rep = os.path.dirname(f)
    for r in map(json.loads, open(f)):
        tag = r["variant"] if r["B"] == 10 else "{}@B{}".format(r["variant"], r["B"])
        row = {"src": "vote", **r, "aln_md5": md5(os.path.join(rep, "vote", tag, "out.fasta"))}
        if "masked_cols" in r:
            row["masked_md5"] = md5(os.path.join(rep, "vote", tag, "out.masked.fasta"))
        rows.append(row)
for r in rows:
    r["identical_to"] = ["{}:{}".format(o["src"], o["variant"]) for o in rows
                         if o is not r and o["rep"] == r["rep"] and o["aln_md5"] and o["aln_md5"] == r["aln_md5"]]
with open(out, "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
