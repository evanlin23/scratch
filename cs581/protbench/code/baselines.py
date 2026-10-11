"""MAFFT L-INS-i alone and Clustal Omega alone (4 threads, 30-min cap) on bbtool_bench work dirs.

    PYTHONPATH=cs581/code python baselines.py WORK OUT.jsonl NAME [NAME ...]

Uses WORK/<NAME>_d0/{unaligned,true}.fasta written by bbtool_bench; scores with
bbtool_bench.acc_ref (estimate restricted to the reference's sequences). A timeout is recorded as such.
"""
import json, os, subprocess, sys, time
from gcmx.bbtool_bench import acc_ref

CAP = int(os.environ.get("CAP", 1800))
TOOLS = {
    "linsi": ["mafft", "--localpair", "--maxiterate", "1000", "--quiet", "--thread", "4", "--anysymbol"],
    "clustalo": ["clustalo", "--threads", "4", "--force", "--outfmt", "fa", "-i"],
}
work, out, names = sys.argv[1], sys.argv[2], sys.argv[3:]
only = os.environ.get("TOOLS", "linsi,clustalo").split(",")
done = set()
if os.path.exists(out):
    done = {(r["dataset"], r["method"]) for r in map(json.loads, open(out))}
for name in names:
    w = os.path.join(work, name + "_d0")
    for tool, argv in TOOLS.items():
        if tool not in only or (name, tool) in done or not os.path.exists(os.path.join(w, "unaligned.fasta")):
            continue
        est = os.path.join(w, "base-" + tool + ".fasta")
        start = time.time()
        try:
            with open(est, "w") as f, open(est + ".log", "w") as e:
                subprocess.run(argv + [os.path.join(w, "unaligned.fasta")], stdout=f, stderr=e, check=True,
                               timeout=CAP)
            row = {"dataset": name, "method": tool, "wall": round(time.time() - start, 1),
                   **acc_ref(os.path.join(w, "true.fasta"), est)}
        except subprocess.TimeoutExpired:
            row = {"dataset": name, "method": tool, "wall": None, "timeout": CAP}
        with open(out, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)
