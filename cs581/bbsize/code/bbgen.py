# AI-assisted (Claude), exploration code for CS581 project
"""New L-INS-i backbones on MAGUS's fixed 25 subsets (MAGUS's sampling: int(size/25) taxa per subset).

    python3 bbgen.py REP OUTDIR SIZE FIRST LAST SEEDBASE [--jobs 4]

REP holds inputs/subalignments; sequences = degapped subset rows (what MAGUS fed its backbones).
Backbone i (FIRST..LAST) uses random.Random(SEEDBASE + i): for each subset in MAGUS's order
(subset_1..subset_25), shuffle its taxa and take the first int(SIZE/25) (all if the subset is smaller).
MAFFT flags = MAGUS's runMafft (--localpair --maxiterate 1000 --ep 0.123 --anysymbol), but --thread 1 and
JOBS backbones at once. Restartable: skips existing OUTDIR/backbone_i_mafft.txt. Timing per backbone in
OUTDIR/timing.jsonl (wall, child CPU = user+sys from wait4).
"""
import argparse, json, os, random, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor


def read_fasta(p):
    d, k = {}, None
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            k = line[1:].strip(); d[k] = []
        elif k is not None:
            d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}


def subsets(rep):
    sd = os.path.join(rep, "inputs", "subalignments")
    fs = sorted(os.listdir(sd), key=lambda f: int(f.split("_")[-1].split(".")[0]))
    out = []
    for f in fs:
        a = read_fasta(os.path.join(sd, f))
        out.append({t: s.replace("-", "").replace(".", "") for t, s in a.items()})
    return out


def run_one(args):
    i, seqs, outdir, size, seed = args
    dst = os.path.join(outdir, "backbone_{}_mafft.txt".format(i))
    if os.path.exists(dst):
        return None
    un = os.path.join(outdir, "backbone_{}_unalign.txt".format(i))
    with open(un, "w") as f:
        for t, s in seqs:
            f.write(">{}\n{}\n".format(t, s))
    tmp = dst + ".tmp"
    t0 = time.time()
    with open(tmp, "w") as out:
        p = subprocess.Popen(["mafft", "--localpair", "--maxiterate", "1000", "--ep", "0.123", "--quiet",
                              "--thread", "1", "--anysymbol", un], stdout=out, stderr=subprocess.DEVNULL)
        _, status, ru = os.wait4(p.pid, 0)
    wall = time.time() - t0
    if status != 0 or os.path.getsize(tmp) == 0:
        raise RuntimeError("mafft failed " + dst)
    os.replace(tmp, dst)
    row = {"i": i, "size": size, "nseq": len(seqs), "seed": seed, "wall": round(wall, 1),
           "cpu": round(ru.ru_utime + ru.ru_stime, 1), "maxrss_mb": round(ru.ru_maxrss / 1024), "jobs": JOBS,
           "load1": round(os.getloadavg()[0], 2)}
    with open(os.path.join(outdir, "timing.jsonl"), "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
    return row


def main():
    global JOBS
    ap = argparse.ArgumentParser()
    ap.add_argument("rep"); ap.add_argument("outdir"); ap.add_argument("size", type=int)
    ap.add_argument("first", type=int); ap.add_argument("last", type=int); ap.add_argument("seedbase", type=int)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    JOBS = a.jobs
    os.makedirs(a.outdir, exist_ok=True)
    subs = subsets(a.rep)
    per = max(1, int(a.size / len(subs)))
    jobs = []
    for i in range(a.first, a.last + 1):
        rng = random.Random(a.seedbase + i)
        seqs = []
        for s in subs:
            taxa = list(s)
            rng.shuffle(taxa)
            seqs += [(t, s[t]) for t in taxa[:per]]
        jobs.append((i, seqs, a.outdir, a.size, a.seedbase + i))
    with ThreadPoolExecutor(a.jobs) as ex:
        list(ex.map(run_one, jobs))


JOBS = 4
if __name__ == "__main__":
    main()
