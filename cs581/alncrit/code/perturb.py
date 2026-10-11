"""Controlled perturbations of a true alignment that raise only SPFN or only SPFP.

    python perturb.py TRUE.fa OUT.fa MODE RATE [--seed 1]

MODE split  (under-alignment): each column with >= 4 residues is, with probability RATE, split
            into two columns by a random half/half partition of its residues. Lost true pairs ->
            SPFN > 0, SPFP = 0, alignment gets longer.
MODE merge  (over-alignment / compression): scanning left to right, with probability RATE an
            adjacent column pair with no row occupied in both is merged into one column. New
            false pairs -> SPFP > 0, SPFN = 0, alignment gets shorter. (On ROSE data the
            compatible columns are sparse insertion columns, so SPFP stays tiny.)
MODE shift  (misplacement): every residue adjacent to a gap in its row is, with probability RATE,
            swapped with that gap (order preserved). Raises SPFN and SPFP together, same length.
"""
import argparse
import random


def read(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return names, ["".join(s) for s in seqs]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("true")
    p.add_argument("out")
    p.add_argument("mode", choices=("split", "merge", "shift"))
    p.add_argument("rate", type=float)
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()
    rng = random.Random(a.seed)
    names, seqs = read(a.true)
    n, L = len(seqs), len(seqs[0])
    cols = ["".join(s[j] for s in seqs) for j in range(L)]
    out = []
    if a.mode == "split":
        for c in cols:
            occ = [i for i, x in enumerate(c) if x != "-"]
            if len(occ) >= 4 and rng.random() < a.rate:
                rng.shuffle(occ)
                half = set(occ[: len(occ) // 2])
                out.append("".join(x if i in half else "-" for i, x in enumerate(c)))
                out.append("".join(x if (i not in half and x != "-") else "-" for i, x in enumerate(c)))
            else:
                out.append(c)
    elif a.mode == "shift":
        rows = [list(x) for x in seqs]
        for r in rows:
            j = 0
            while j < L - 1:
                if (r[j] != "-") != (r[j + 1] != "-") and rng.random() < a.rate:
                    r[j], r[j + 1] = r[j + 1], r[j]
                    j += 2
                else:
                    j += 1
        out = ["".join(r[j] for r in rows) for j in range(L)]
    else:
        cur = cols[0]
        for c in cols[1:]:
            compatible = all(x == "-" or y == "-" for x, y in zip(cur, c))
            if compatible and any(x != "-" for x in cur) and any(y != "-" for y in c) and rng.random() < a.rate:
                cur = "".join(x if x != "-" else y for x, y in zip(cur, c))
            else:
                out.append(cur)
                cur = c
        out.append(cur)
    with open(a.out, "w") as f:
        for i, name in enumerate(names):
            f.write(">%s\n%s\n" % (name, "".join(c[i] for c in out)))


if __name__ == "__main__":
    main()
