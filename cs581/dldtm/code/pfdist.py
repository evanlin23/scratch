"""Phyloformer (Nesterenko et al., MBE 2025, doi:10.1093/molbev/msaf051) as a CPU
distance estimator for DNA alignments, followed by FastME (BME + NNI + SPR).

The public checkpoints are protein models (22-letter alphabet ARNDCQEGHILKMFPSTWYVX-).
DNA is fed in by mapping A,C,G,T to the amino-acid channels with the same letters
(A=Ala, C=Cys, G=Gly, T=Thr); anything else becomes X and '-' stays a gap. With the
pretrained weights this is out of domain; `train_nt.py` fine-tunes on DNA encoded the
same way, so the same code serves both.

Long alignments are cut into windows of `win` sites (the model was trained on ~200-500
sites). Distances are averaged over windows, weighted by window length.

Usage as a script: python3 pfdist.py CKPT ALN.fasta OUT.tre [win]
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import torch

PF_REPO = os.environ.get("PF_REPO", "/opt/src/Phyloformer")
sys.path.insert(0, PF_REPO)
from phyloformer.model import Phyloformer  # noqa: E402

FASTME = os.environ.get("FASTME", f"{PF_REPO}/bin/bin_linux/fastme")
ALPHABET = "ARNDCQEGHILKMFPSTWYVX-"
LUT = np.full(256, ALPHABET.index("X"), dtype=np.int64)
for c in ALPHABET:
    LUT[ord(c)] = ALPHABET.index(c)
    LUT[ord(c.lower())] = ALPHABET.index(c)
# DNA: only A/C/G/T/gap are meaningful; everything else (N, IUPAC, U?) -> X
for c in "RNDQEHILKMFPSWYVB":
    LUT[ord(c)] = LUT[ord(c.lower())] = ALPHABET.index("X")
LUT[ord("U")] = LUT[ord("u")] = ALPHABET.index("T")
LUT[ord(".")] = ALPHABET.index("-")

_cache = {}


def load_model(ckpt):
    if ckpt in _cache:
        return _cache[ckpt]
    c = torch.load(ckpt, map_location="cpu", weights_only=False)
    params = dict(c["hyper_parameters"])
    params["device"] = "cpu"
    m = Phyloformer(**params)
    sd = c["state_dict"] if "state_dict" in c else c
    m.load_state_dict({k.replace("model.", ""): v for k, v in sd.items() if k != "model.seq2pair"}, strict=False)
    m.eval()
    _cache[ckpt] = m
    return m


def encode(seqs):
    """seqs: list of equal-length strings -> LongTensor (L, n) of alphabet indices."""
    a = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), -1)
    return torch.from_numpy(LUT[a].T.copy())


def onehot(idx):
    """(L, n) indices -> (1, 22, L, n) float."""
    return torch.nn.functional.one_hot(idx, num_classes=22).permute(2, 0, 1).unsqueeze(0).float()


def distances(model, names, seqs, win=256):
    """Return n x n numpy distance matrix (expected substitutions per site)."""
    idx = encode(seqs)
    # drop all-gap columns of this subset
    keep = (idx != ALPHABET.index("-")).any(dim=1)
    idx = idx[keep]
    L, n = idx.shape
    nwin = max(1, round(L / win))
    bounds = np.linspace(0, L, nwin + 1).astype(int)
    acc = None
    tot = 0
    with torch.no_grad():
        for a, b in zip(bounds[:-1], bounds[1:]):
            out = model(onehot(idx[a:b])).reshape(-1)
            acc = out * (b - a) if acc is None else acc + out * (b - a)
            tot += b - a
    v = (acc / tot).numpy()
    D = np.zeros((n, n))
    iu = np.triu_indices(n, 1)
    D[iu] = v
    return D + D.T


def write_phylip(D, names, path):
    with open(path, "w") as f:
        f.write(f"{len(names)}\n")
        for i, nm in enumerate(names):
            f.write(nm + " " + " ".join("%.8f" % x for x in D[i]) + "\n")


def fastme_tree(D, names, out):
    """BME tree with NNI+SPR (Phyloformer paper's protocol). Names are mapped to
    short ids because FastME truncates long labels."""
    with tempfile.TemporaryDirectory() as td:
        ids = ["s%d" % i for i in range(len(names))]
        write_phylip(D, ids, f"{td}/d.phy")
        subprocess.run([FASTME, "-i", f"{td}/d.phy", "-o", f"{td}/t.nwk", "--nni", "--spr", "-T", "1"],
                       check=True, capture_output=True)
        s = open(f"{td}/t.nwk").read().strip()
    import re
    s = re.sub(r"\bs(\d+)\b(?=[:,)])", lambda m: names[int(m.group(1))], s)
    with open(out, "w") as f:
        f.write(s + "\n")


def read_fasta(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            seqs.append([])
        elif line:
            seqs[-1].append(line)
    return names, ["".join(s).upper() for s in seqs]


if __name__ == "__main__":
    ckpt, aln, out = sys.argv[1:4]
    win = int(sys.argv[4]) if len(sys.argv) > 4 else 256
    names, seqs = read_fasta(aln)
    D = distances(load_model(ckpt), names, seqs, win)
    fastme_tree(D, names, out)
