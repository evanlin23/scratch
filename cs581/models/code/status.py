"""Print one line per (dataset, method) from worker logs."""
import json, sys
for p in sys.argv[1:]:
    for l in open(p):
        try:
            d = json.loads(l)
        except ValueError:
            print(l.strip()[:150]); continue
        for k, v in d.items():
            for m, r in v.items():
                print(k, m, {x: (round(y, 3) if isinstance(y, float) else y) for x, y in r.items() if x in ("wall", "avgErr", "treeFN", "error")})
