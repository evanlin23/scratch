"""Drop non-ok rows of TOOL from a results jsonl (so an interrupted step can be rerun).

    python3 drop_rows.py RESULTS.jsonl TOOL
"""

import json
import sys

path, tool = sys.argv[1], sys.argv[2]
rows = [l for l in open(path) if not (json.loads(l)["tool"] == tool and json.loads(l)["status"] != "ok")]
open(path, "w").writelines(rows)
