"""Run the method atlas on disjoint blocks of a 4-taxon adversarial pool (e.g. the
ASTRAL-Pro counterexample cand2): fraction of datasets with a wrong species tree.
usage: python atlas4.py POOL.nwk SPECIES_NWK OUT.jsonl KS [methods]"""
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core  # noqa: E402
from phylo import parse_newick  # noqa: E402
from run_curve import ATLAS, method_tree  # noqa: E402

pool, sp_nwk, out, ks = sys.argv[1], sys.argv[2], sys.argv[3], [int(x) for x in sys.argv[4].split(",")]
methods = sys.argv[5].split(",") if len(sys.argv) > 5 else ATLAS
fams = [l.strip() for l in open(pool) if ";" in l]
st = parse_newick(sp_nwk)
species = sorted(st.label[v] for v in st.leaves())
spi = {s: i for i, s in enumerate(species)}
done = set()
if os.path.exists(out):
    done = {(r["method"], r["K"], r["block"]) for r in map(json.loads, open(out))}
os.makedirs("/opt/runs/gdlcons/tmp", exist_ok=True)
for k in ks:
    for b in range(min(len(fams) // k, 20)):
        todo = [m for m in methods if (m, k, b) not in done]
        if not todo:
            continue
        rts = [core.true_rt(f, spi)[0] for f in fams[b * k:(b + 1) * k]]
        plain = [rt.newick() for rt in rts]
        with tempfile.TemporaryDirectory(dir="/opt/runs/gdlcons/tmp") as td:
            for m in todo:
                rec = {"method": m, "K": k, "block": b}
                try:
                    est = method_tree(m, plain, rts, species, td)
                    rec["FN"] = core.rf_error(est, st)[0]
                except Exception as e:  # noqa: BLE001
                    rec["error"] = "%s: %s" % (type(e).__name__, str(e)[:200])
                with open(out, "a") as f:
                    f.write(json.dumps(rec) + "\n")
    print("K", k, "done", flush=True)
