"""Minimal FASTA helpers (ordered dicts of name -> sequence)."""


def read(path):
    seqs, name, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def write(seqs, path, width=0):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(">{}\n".format(name))
            if width:
                for i in range(0, len(seq), width):
                    f.write(seq[i:i + width] + "\n")
            else:
                f.write(seq + "\n")


def ungap(seqs):
    return {n: s.replace("-", "").replace(".", "") for n, s in seqs.items()}


def restrict(aln, taxa):
    """Induced sub-alignment on `taxa`, with all-gap columns removed."""
    rows = [aln[t] for t in taxa]
    keep = [i for i in range(len(rows[0])) if any(r[i] not in "-." for r in rows)]
    return {t: "".join(aln[t][i] for i in keep) for t in taxa}


def upper(seqs):
    return {n: s.upper() for n, s in seqs.items()}
