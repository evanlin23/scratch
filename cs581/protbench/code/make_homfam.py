"""Subsample a HomFam family to N sequences, always keeping the Homstrad seed (reference) sequences.

    python make_homfam.py FAMILY_DIR OUT_DIR N SEED

FAMILY_DIR holds all.unaln.fasta (whole family) and all.aln.fasta (seed reference), as in the SALMA/EMMA
data release (Illinois Data Bank IDB-2567453, homfam/<family>/<replicate>/). Writes OUT_DIR/unaln.fasta
(N sequences: the seeds plus a uniform random sample of the rest) and OUT_DIR/ref.fasta (the seeds).
"""
import os, random, sys


def read(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write(seqs, path):
    with open(path, "w") as f:
        for k, v in seqs.items():
            f.write(">{}\n{}\n".format(k, v))


fam, out, n, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
allu = read(os.path.join(fam, "all.unaln.fasta"))
ref = read(os.path.join(fam, "all.aln.fasta"))
for k, v in ref.items():
    u = v.replace("-", "").replace(".", "").upper()
    assert k in allu and allu[k].upper() == u, k
rest = sorted(k for k in allu if k not in ref)
random.Random(seed).shuffle(rest)
keep = list(ref) + rest[:n - len(ref)]
os.makedirs(out, exist_ok=True)
write({k: allu[k] for k in keep}, os.path.join(out, "unaln.fasta"))
write(ref, os.path.join(out, "ref.fasta"))
print(fam, len(allu), "->", len(keep), "seeds", len(ref))
