"""Pretty-print score.json lines: python show.py DS [DS ...]"""
import json, os, sys
for ds in sys.argv[1:]:
    for m in sorted(os.listdir(ds)):
        p = os.path.join(ds, m, "score.json")
        if os.path.exists(p):
            d = json.load(open(p))
            print(os.path.basename(ds), m, round(d["time"]), " ".join("%s=%.3f/%.3f" % (k, v["SPFN"], v["SPFP"]) for k, v in d["score"].items()))
