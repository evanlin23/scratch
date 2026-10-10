#!/usr/bin/env python3
"""TIPP3-fast style BLAST alignment: map each query through the HSP of its top BLASTN hit
(default task = megablast, as TIPP3 calls blastn) onto that backbone sequence's columns.

usage: blast_aln.py INST OUTDIR THREADS [TASK]   (TASK: megablast = TIPP3 default, or blastn = sensitive)
writes OUTDIR/blast_map.json (query -> per-residue backbone column or -1), OUTDIR/tophit.json
(query -> (backbone name, bitscore)) and OUTDIR/blast_time.txt
"""
import sys, json, time, os, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score import read_fasta

BIN = '/opt/mm/root/envs/blast/bin'


def rc(s):
    m = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(m.get(c, c) for c in reversed(s))


def main():
    inst, out, th = sys.argv[1], sys.argv[2], sys.argv[3]
    task = sys.argv[4] if len(sys.argv) > 4 else 'megablast'
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    db = f'{out}/bbdb'
    subprocess.run([f'{BIN}/makeblastdb', '-in', f'{inst}/backbone.unaln.fasta', '-dbtype', 'nucl',
                    '-out', db], check=True, stdout=subprocess.DEVNULL)
    t1 = time.time()
    res = subprocess.run([f'{BIN}/blastn', '-db', db, '-query', f'{inst}/queries.fasta',
                          '-outfmt', '6 qseqid sseqid qstart qend sstart send qseq sseq bitscore',
                          '-max_target_seqs', '5', '-num_threads', th, '-task', task],
                         check=True, capture_output=True, text=True).stdout
    t2 = time.time()
    bb = read_fasta(f'{inst}/backbone.fasta')
    qs = read_fasta(f'{inst}/queries.fasta')
    colidx = {}
    for n, s in bb.items():
        colidx[n] = [c for c, ch in enumerate(s) if ch != '-']
    maps, top, top5 = {}, {}, {}
    for line in res.splitlines():
        q, sid, qs_, qe, ss, se, qa, sa, bits = line.split('\t')
        hits = top5.setdefault(q, [])
        if sid not in [h[0] for h in hits]:
            hits.append((sid, float(bits)))
        if q in maps:
            continue  # keep the first (top) HSP only
        qs_, qe, ss, se = int(qs_), int(qe), int(ss), int(se)
        L = len(qs[q])
        m = [-1] * L
        cols = colidx[sid]
        if ss <= se:
            qi, ti = qs_ - 1, ss - 1
            for a, b in zip(qa, sa):
                if a != '-' and b != '-':
                    m[qi] = cols[ti]
                if a != '-':
                    qi += 1
                if b != '-':
                    ti += 1
        else:
            # minus strand: query positions increase while subject positions decrease;
            # queries are all forward-strand here, so this branch only maps what it can
            qi, ti = qs_ - 1, ss - 1
            for a, b in zip(qa, sa):
                if a != '-' and b != '-':
                    m[qi] = -1
                if a != '-':
                    qi += 1
                if b != '-':
                    ti -= 1
        maps[q] = m
        top[q] = (sid, float(bits))
    for q in qs:
        maps.setdefault(q, [-1] * len(qs[q]))
    t3 = time.time()
    json.dump(maps, open(f'{out}/blast_map.json', 'w'))
    json.dump(top, open(f'{out}/tophit.json', 'w'))
    json.dump(top5, open(f'{out}/tophits5.json', 'w'))
    with open(f'{out}/blast_time.txt', 'w') as f:
        f.write(f'makeblastdb {t1-t0:.2f}\nblastn {t2-t1:.2f}\nextract {t3-t2:.2f}\ntotal {t3-t0:.2f}\n')
    print(f'blast total {t3-t0:.1f}s, {len(top)}/{len(qs)} queries with hits')


if __name__ == '__main__':
    main()
