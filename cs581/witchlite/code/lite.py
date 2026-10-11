#!/usr/bin/env python3
"""WITCH-lite HMM selection: score each query against a subset of the WITCH/UPP HMM ensemble,
compute WITCH weights (adjusted bit-scores) over the scored HMMs only, and write a WITCH
weights.txt so that WITCH's own query-alignment + merging stage can run unchanged.

usage: lite.py INST HMMDIR OUTDIR STRATEGY [--k 10] [--tau 1.0] [--beam 2] [--threads 4]
  STRATEGY: all | hier | hier_es | beam | blastpath | blastpath_sib
  --tau t < 1: adaptive k (keep the fewest top-weighted HMMs whose weights sum to >= t, at most k)
writes OUTDIR/weights.txt, OUTDIR/scoring.json (time, #HMM-query scores)
"""
import os, sys, re, json, time, glob, argparse
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import subprocess, threading, tempfile


def load_ensemble(hmmdir):
    dirs = [d for d in glob.glob(f'{hmmdir}/**/A_0_*', recursive=True) if os.path.isdir(d)]
    nodes = {}
    for d in dirs:
        idx = int(os.path.basename(d).split('_')[2])
        parent = None
        model = [f for f in glob.glob(f'{d}/hmmbuild.model.*')]
        inp = [f for f in glob.glob(f'{d}/hmmbuild.input.*.fasta')]
        size = sum(1 for l in open(inp[0]) if l.startswith('>'))
        names = [l[1:].split()[0] for l in open(inp[0]) if l.startswith('>')]
        nodes[idx] = dict(dir=d, parent=parent, size=size, model=model[0], names=set(names), children=[])
    # WITCH writes the (nested) decomposition flat on disk: recover parent = smallest strict superset
    order = sorted(nodes, key=lambda i: nodes[i]['size'])
    for a_i, i in enumerate(order):
        for j in order[a_i + 1:]:
            if nodes[j]['size'] > nodes[i]['size'] and nodes[i]['names'] <= nodes[j]['names']:
                nodes[i]['parent'] = j
                break
    for i, n in nodes.items():
        if n['parent'] is not None:
            nodes[n['parent']]['children'].append(i)
    roots = [i for i, n in nodes.items() if n['parent'] is None]
    return nodes, roots


class Scorer:
    def __init__(self, nodes, seqs, threads, tmp):
        self.nodes, self.seqs, self.threads, self.tmp = nodes, seqs, threads, tmp
        self.n_scores = 0

    def _one(self, args):
        # same binary and flags as WITCH's all-against-all search (hmmsearch --max -E 99999999),
        # restricted to the queries routed to this HMM
        idx, qids = args
        fa = f'{self.tmp}/q.{idx}.{threading.get_ident()}.fa'
        tb = fa + '.tbl'
        with open(fa, 'w') as f:
            for q in qids:
                f.write(f'>{q}\n{self.seqs[q]}\n')
        subprocess.run(['hmmsearch', '--cpu', '1', '--noali', '-E', '99999999', '--max', '-o', '/dev/null',
                        '--tblout', tb, self.nodes[idx]['model'], fa], check=True)
        out = {}
        for line in open(tb):
            if line.startswith('#'):
                continue
            f = line.split()
            out[f[0]] = float(f[5])
        os.remove(fa); os.remove(tb)
        return idx, qids, out

    def score(self, jobs):
        """jobs: dict hmm_idx -> list of query names. returns dict (q, idx) -> bitscore."""
        res = {}
        # split big batches (e.g. the root HMM, scored for every query) so all threads are busy,
        # as WITCH splits its queries into num_cpus chunks
        jl = []
        for i, q in jobs.items():
            step = max(50, -(-len(q) // self.threads))
            jl.extend((i, q[j:j + step]) for j in range(0, len(q), step))
        # largest jobs first for load balance
        jl.sort(key=lambda x: -len(x[1]))
        with ThreadPoolExecutor(self.threads) as ex:
            for idx, qids, out in ex.map(self._one, jl):
                self.n_scores += len(qids)
                for q in qids:
                    res[(q, idx)] = out.get(q, -1e9)
        return res


def adj(bits, size):
    return bits + np.log2(size)


def weights_from_scores(scored, nodes, k, tau):
    """scored: dict idx -> bitscore for one query. WITCH weight w_i = s_i 2^b_i / sum_j s_j 2^b_j."""
    idx = list(scored.keys())
    a = np.array([adj(scored[i], nodes[i]['size']) for i in idx])
    a = a - a.max()
    w = np.power(2.0, a)
    w = w / w.sum()
    order = np.argsort(-w)[:k]
    out = []
    cum = 0.0
    for o in order:
        out.append((idx[o], float(w[o])))
        cum += w[o]
        if cum >= tau:
            break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('inst'); ap.add_argument('hmmdir'); ap.add_argument('outdir'); ap.add_argument('strategy')
    ap.add_argument('--k', type=int, default=10); ap.add_argument('--tau', type=float, default=1.0)
    ap.add_argument('--beam', type=int, default=2); ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--blastdir', default=None)
    ap.add_argument('--cache', default=None, help='reuse scores from a previous "all" run (selection-only experiments)')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    t_load0 = time.time()
    nodes, roots = load_ensemble(a.hmmdir)
    seqs = {}
    for line in open(f'{a.inst}/queries.fasta'):
        line = line.strip()
        if line.startswith('>'):
            nm = line[1:].split()[0]; seqs[nm] = ''
        elif line:
            seqs[nm] += line
    qnames = list(seqs.keys())
    t_load = time.time() - t_load0

    sc = Scorer(nodes, seqs, a.threads, tempfile.mkdtemp(dir=a.outdir))
    t0 = time.time()
    per_q = {q: {} for q in qnames}

    def run(jobs):
        r = sc.score(jobs)
        for (q, i), b in r.items():
            per_q[q][i] = b

    if a.strategy == 'all':
        run({i: qnames for i in nodes})
    elif a.strategy in ('hier', 'hier_es', 'beam'):
        beam = a.beam if a.strategy == 'beam' else 1
        run({r: qnames for r in roots})
        frontier = {q: sorted(roots, key=lambda r: -adj(per_q[q][r], nodes[r]['size']))[:beam] for q in qnames}
        while True:
            jobs = {}
            for q, fr in frontier.items():
                for f in fr:
                    for c in nodes[f]['children']:
                        jobs.setdefault(c, []).append(q)
            if not jobs:
                break
            run(jobs)
            newf = {}
            for q, fr in frontier.items():
                ch = [c for f in fr for c in nodes[f]['children']]
                if not ch:
                    continue
                ch.sort(key=lambda c: -adj(per_q[q][c], nodes[c]['size']))
                if a.strategy == 'hier_es':
                    cur = adj(per_q[q][fr[0]], nodes[fr[0]]['size'])
                    if adj(per_q[q][ch[0]], nodes[ch[0]]['size']) < cur:
                        continue
                newf[q] = ch[:beam]
            frontier = newf
    elif a.strategy in ('blastpath', 'blastpath_sib'):
        top = json.load(open(f'{a.blastdir}/tophit.json'))
        # deepest node containing each backbone sequence -> path to root
        leaf_of = {}
        depth = {}
        def dep(i):
            if i not in depth:
                p = nodes[i]['parent']
                depth[i] = 0 if p is None else dep(p) + 1
            return depth[i]
        for i, n in nodes.items():
            for nm in n['names']:
                if nm not in leaf_of or dep(i) > dep(leaf_of[nm]):
                    leaf_of[nm] = i
        jobs = {}
        for q in qnames:
            if q in top:
                path = []
                i = leaf_of[top[q][0]]
                while i is not None:
                    path.append(i); i = nodes[i]['parent']
                sel = set(path)
                if a.strategy == 'blastpath_sib':
                    for i in path:
                        p = nodes[i]['parent']
                        if p is not None:
                            sel.update(nodes[p]['children'])
            else:
                sel = set(nodes)  # no BLAST hit: fall back to scoring everything
            for i in sel:
                jobs.setdefault(i, []).append(q)
        run(jobs)
    else:
        raise SystemExit('unknown strategy')
    t_score = time.time() - t0

    # weights for several (k, tau) settings from the same scores (scoring time is shared)
    settings = [(a.k, a.tau), (10, 1.0), (10, 0.99), (10, 0.95), (3, 1.0), (1, 1.0)]
    mean_k = {}
    for k, tau in dict.fromkeys(settings):
        tag = f'k{k}' + ('' if tau >= 1 else f'_t{int(round(tau * 100))}')
        ks = []
        with open(f'{a.outdir}/weights_{tag}.txt', 'w') as f:
            for q in qnames:
                w = weights_from_scores(per_q[q], nodes, k, tau)
                ks.append(len(w))
                f.write('{}:{}\n'.format(q, tuple(w)))
        mean_k[tag] = float(np.mean(ks))
    info = dict(strategy=a.strategy, beam=a.beam, n_hmms=len(nodes),
                n_queries=len(qnames), n_scores=sc.n_scores, scores_per_query=sc.n_scores / len(qnames),
                score_time=t_score, load_time=t_load, mean_k=mean_k)
    if a.strategy == 'all':
        json.dump({q: {str(i): b for i, b in d.items()} for q, d in per_q.items()},
                  open(f'{a.outdir}/allscores.json', 'w'))
    json.dump(info, open(f'{a.outdir}/scoring.json', 'w'), indent=1)
    print(json.dumps(info))


if __name__ == '__main__':
    main()
