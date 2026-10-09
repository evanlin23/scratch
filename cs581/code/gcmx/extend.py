"""HMM-extend backbone alignments to all sequences (MAGUS "Slow" mode, precomputed).

    python -m gcmx.extend BACKBONE_DIR UNALIGNED.fa OUT_DIR [--jobs 4]

For each backbone: hmmbuild on it, hmmalign every other sequence to the HMM,
and write backbone + extended sequences as one FASTA. Extended sequences keep
HMMER's insertion convention (lower-case letters, '.'), which MAGUS's graph
builder already understands (insertions are not used as evidence). Uses the
same binaries and flags as MAGUS's --graphbuildhmmextend true, but runs once
per replicate so that different merges can share the same extended evidence.
"""

import argparse
import concurrent.futures
import os
import subprocess
import tempfile

from magus.helpers import sequenceutils

from . import fasta

TOOLS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "MAGUS", "magus", "tools", "hmmer")


def extend(backbone_path, unaligned, out_path):
    backbone = fasta.read(backbone_path)
    queries = {n: s for n, s in unaligned.items() if n not in backbone}
    with tempfile.TemporaryDirectory() as tmp:
        hmm = os.path.join(tmp, "model.hmm")
        subprocess.run([os.path.join(TOOLS, "hmmbuild"), "--ere", "0.59", "--cpu", "1", "--symfrac", "0.0",
                        "--informat", "afa", hmm, backbone_path], check=True, capture_output=True)
        query_path = os.path.join(tmp, "queries.fa")
        fasta.write(queries, query_path)
        sto = os.path.join(tmp, "aligned.sto")
        subprocess.run([os.path.join(TOOLS, "hmmalign"), "-o", sto, hmm, query_path], check=True, capture_output=True)
        extended = sequenceutils.readFromStockholm(sto, includeInsertions=True)
    with open(out_path, "w") as f:
        for name, seq in backbone.items():
            f.write(">{}\n{}\n".format(name, seq))
        for name, record in extended.items():
            f.write(">{}\n{}\n".format(name, record.seq))
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("backbones")
    parser.add_argument("unaligned")
    parser.add_argument("outdir")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    unaligned = fasta.upper(fasta.read(args.unaligned))
    os.makedirs(args.outdir, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(extend, os.path.join(args.backbones, f), unaligned, os.path.join(args.outdir, f))
                   for f in sorted(os.listdir(args.backbones))]
        for fut in futures:
            print("extended", fut.result())


if __name__ == "__main__":
    main()
