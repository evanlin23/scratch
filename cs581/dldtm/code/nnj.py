"""NeuralNJ subset-tree wrapper: FASTA -> PHYLIP with names t1..tn -> nnj_driver.py -> Newick."""
import os
import re
import subprocess
import tempfile

import common

PY = os.environ.get("NNJ_PY", "/opt/mm/root/envs/nnj/bin/python")
HERE = os.path.dirname(os.path.abspath(__file__))


def tree(aln, out):
    s = common.read_fasta(aln)
    names = list(s)
    with tempfile.TemporaryDirectory() as td:
        with open(f"{td}/a.phy", "w") as f:
            f.write(f"{len(s)} {len(s[names[0]])}\n")
            for i, n in enumerate(names):
                f.write(f"t{i + 1} {s[n]}\n")
        subprocess.run([PY, f"{HERE}/nnj_driver.py", f"{td}/a.phy", f"{td}/t.nwk"], check=True,
                       capture_output=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
        nw = open(f"{td}/t.nwk").read().strip()
    nw = re.sub(r"\bt(\d+)\b(?=[:,);])", lambda m: names[int(m.group(1)) - 1], nw)
    with open(out, "w") as f:
        f.write(nw + "\n")
