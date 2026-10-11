"""Full-data IQ-TREE (-fast, GTR+G) outside pipe.py's lock (for running beside it).
Usage: python3 iq_full.py DIR  -> DIR/full_iq.json  (env THREADS, IQ_CAP seconds, IQ_MEM)"""
import json, os, sys
from pipe import timed, IQ, THREADS, IQ_CAP, read_tree, fn_fp
d = sys.argv[1]
w, m, rc = timed([IQ, "-s", f"{d}/aln.fa", "-m", "GTR+G", "-T", str(THREADS), "-fast", "-seed", "1", "--prefix",
                  f"{d}/full_iq", "-redo", "-quiet", "-mem", os.environ.get("IQ_MEM", "10G")], timeout=IQ_CAP)
r = dict(wall=w, mem=m, rc=rc)
if os.path.exists(f"{d}/full_iq.treefile"):
    r["fn"] = fn_fp(read_tree(f"{d}/full_iq.treefile"), read_tree(f"{d}/true.tre"))[0]
json.dump(r, open(f"{d}/full_iq.json", "w"), indent=1)
print(r)
