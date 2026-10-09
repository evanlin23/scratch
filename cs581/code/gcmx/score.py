"""Score an estimated alignment against a reference with FastSP.

    python -m gcmx.score REFERENCE ESTIMATE  ->  prints a JSON dict

Both alignments are upper-cased first (FastSP can treat lower case
specially). Reported: SPFN, SPFP, expansion/compression, TC.
"""

import json
import os
import subprocess
import sys
import tempfile

from . import fasta

FASTSP_JAR = os.environ.get("FASTSP_JAR", "/opt/tools/FastSP/FastSP.jar")


def fastsp(reference, estimate):
    with tempfile.TemporaryDirectory() as tmp:
        ref = os.path.join(tmp, "ref.fa")
        est = os.path.join(tmp, "est.fa")
        fasta.write(fasta.upper(fasta.read(reference)), ref)
        fasta.write(fasta.upper(fasta.read(estimate)), est)
        out = subprocess.run(["java", "-Xmx4g", "-jar", FASTSP_JAR, "-r", ref, "-e", est],
                             capture_output=True, text=True, check=True)
    stats = {}
    for line in (out.stdout + out.stderr).splitlines():
        tokens = line.split()
        if len(tokens) == 2 and tokens[0] in ("SP-Score", "Modeler", "SPFN", "SPFP", "Compression", "TC"):
            stats[tokens[0]] = float(tokens[1])
        elif line.startswith("MaxLenNoGap"):
            fields = dict(f.strip().split("= ") for f in line.split(","))
            stats["LenRef"], stats["LenEst"] = int(fields["LenRef"]), int(fields["LenEst"])
    stats["avgErr"] = (stats["SPFN"] + stats["SPFP"]) / 2
    return stats


if __name__ == "__main__":
    print(json.dumps(fastsp(sys.argv[1], sys.argv[2])))
