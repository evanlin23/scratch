#!/usr/bin/env python3
"""Restartable driver: WITCH default, BLAST, lite selectors, WITCH alignment stage per variant, scoring.
usage: run_inst.py INST RESULTS_JSON
Every step is skipped if its output exists, so the script can be relaunched after a kill.
"""
import os, sys, json, subprocess, time, re
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from score import load_instance, map_from_alignment, metrics

inst, resf = sys.argv[1], sys.argv[2]
H = f'{inst}/witch_default/tree_decomp/root'
res = json.load(open(resf)) if os.path.exists(resf) else {}


def save():
    json.dump(res, open(resf, 'w'), indent=1)


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


# 1) WITCH default (decomposition + all-vs-all hmmsearch + alignment)
if not os.path.exists(f'{inst}/witch_default/aln.fasta'):
    s = time.time()
    sh(f'cd {inst} && witch.py -b {inst}/backbone.fasta -e {inst}/backbone.tre -q {inst}/queries.fasta '
       f'-d {inst}/witch_default -o aln.fasta -t 4 --molecule dna --save-weight 1 --keep-decomposition 1 '
       f'> {inst}/witch_default.log 2>&1')
    open(f'{inst}/witch_default.wall', 'w').write(f'wall {time.time() - s}\n')
rb = open(f'{inst}/witch_default/runtime_breakdown.txt').read()
retimed = os.path.exists(f'{inst}/retime.wall')
if retimed:  # clean re-run of WITCH default on the idle machine (the first run overlapped another job)
    rb = open(f'{inst}/witch_retime/runtime_breakdown.txt').read()
dec = float(re.search(r'decompose the backbone \(s\): ([\d.]+)', rb).group(1))
srch = float(re.search(r'HMMSearches \(s\): ([\d.]+)', rb).group(1))
wall = float(open(f'{inst}/retime.wall').read()) if retimed else float(open(f'{inst}/witch_default.wall').read().split()[1])
res['witch_default_time'] = dict(wall=wall, decomposition=dec, search=srch, rest=wall - dec - srch)

truth, cc, bb = load_instance(inst)
bbn = list(bb)


def record(name, est, extra):
    fn, fp, fnq, fpq = metrics(truth, cc, est)
    res.setdefault('methods', {})[name] = dict(spfn=fn, spfp=fp, spfn_q=fnq.tolist(), spfp_q=fpq.tolist(), **extra)
    save()
    print(f'{name:28s} SPFN {fn:.4f} SPFP {fp:.4f} {extra}', flush=True)


est, _ = map_from_alignment(f'{inst}/witch_default/aln.fasta', bbn, truth)
record('WITCH', est, dict(time=wall, time_nodecomp=wall - dec))

# 2) BLAST (TIPP3-fast alignment)
if not os.path.exists(f'{inst}/blast/blast_map.json'):
    sh(f'python3 {HERE}/blast_aln.py {inst} {inst}/blast 4')
bt = float(open(f'{inst}/blast/blast_time.txt').read().split('total')[1])
record('BLAST', json.load(open(f'{inst}/blast/blast_map.json')), dict(time=bt, time_nodecomp=bt))
# sensitive BLASTN (word size 11); also the guide for the blastpath selectors
if not os.path.exists(f'{inst}/blast_sens/blast_map.json'):
    sh(f'python3 {HERE}/blast_aln.py {inst} {inst}/blast_sens 4 blastn')
bt = float(open(f'{inst}/blast_sens/blast_time.txt').read().split('total')[1])
record('BLAST-sens', json.load(open(f'{inst}/blast_sens/blast_map.json')), dict(time=bt, time_nodecomp=bt))

# 3) selectors
for strat in ['all', 'hier', 'hier_es', 'beam', 'blastpath', 'blastpath_sib', 'beam3', 'beam4']:
    od = f'{inst}/lite_{strat}'
    if not os.path.exists(f'{od}/scoring.json'):
        extra = f'--blastdir {inst}/blast_sens' if strat.startswith('blastpath') else ''
        if strat.startswith('beam') and strat != 'beam':
            extra = f'--beam {strat[4:]}'; strat = 'beam'
        sh(f'python3 {HERE}/lite.py {inst} {H} {od} {strat} {extra} > /dev/null')
# hybrid rule (BLAST path(s)+siblings if top hit >= 100 bits, else beam-3), replayed from cached scores
if not os.path.exists(f'{inst}/sim/weights_hyb_b3_100.txt'):
    sh(f'python3 {HERE}/sim.py {inst} {inst}/sim --write hyb_b3_100 > {inst}/sim.log')

# 4) WITCH alignment stage on each selection
variants = [('all', 'k10'), ('all', 'k1'), ('all', 'k3'), ('all', 'k10_t99'), ('all', 'k10_t95'),
            ('hier', 'k10'), ('hier', 'k10_t99'), ('hier_es', 'k10'), ('beam', 'k10'), ('beam', 'k10_t99'),
            ('blastpath', 'k10'), ('blastpath', 'k10_t99'), ('blastpath_sib', 'k10'), ('blastpath_sib', 'k10_t99'),
            ('blastpath_sib', 'k10_t95'), ('beam3', 'k10'), ('beam3', 'k10_t99'), ('beam4', 'k10'),
            ('hyb', 'b3_100')]
for strat, tag in variants:
    st = f'{inst}/stage_{strat}_{tag}'
    if strat == 'hyb':
        # selection time ESTIMATED: sensitive BLAST + (HMM scores/query) x (per-score cost of online beam-3)
        sim = json.load(open(f'{inst}/sim/sim.json'))['hyb_b3_100']
        b3 = json.load(open(f'{inst}/lite_beam3/scoring.json'))
        sc = dict(scores_per_query=sim['scores_per_query'], mean_k={tag: 10.0},
                  score_time=sim['scores_per_query'] * b3['score_time'] / b3['scores_per_query'])
        wfile, k = f'{inst}/sim/weights_hyb_b3_100.txt', 10
    else:
        od = f'{inst}/lite_{strat}'
        sc = json.load(open(f'{od}/scoring.json'))
        wfile, k = f'{od}/weights_{tag}.txt', tag.split('_')[0][1:]
    if not os.path.exists(f'{st}/aln.fasta'):
        sh(f'bash {HERE}/witch_stage.sh {inst} {H} {wfile} {st} {k}')
    stage = float(open(f'{st}/stage_wall.txt').read())
    selt = sc['score_time'] + (bt if strat.startswith('blastpath') or strat == 'hyb' else 0)
    est, _ = map_from_alignment(f'{st}/aln.fasta', bbn, truth)
    record(f'{strat}/{tag}', est, dict(select_time=selt, stage_time=stage, time=dec + selt + stage,
                                       time_nodecomp=selt + stage, scores_per_query=sc['scores_per_query'],
                                       mean_k=sc['mean_k'][tag]))
print('done')
