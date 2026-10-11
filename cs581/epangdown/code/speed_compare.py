"""Speed-test placement agreement: share of queries whose best edge matches between EPA-ng builds
(same reference tree, so edge numbers are comparable). Usage: python speed_compare.py <speeddir>"""
import sys, json, glob, os, itertools
def best(p):
    J = json.load(open(p)); f = J['fields']; ie, il = f.index('edge_num'), f.index('like_weight_ratio')
    return {(x.get('n') or [m[0] for m in x['nm']])[0]: max(x['p'], key=lambda r: r[il])[ie] for x in J['placements']}
for kd in sorted(glob.glob(f'{sys.argv[1]}/k*')):
    for q in ('frag', 'full'):
        B = {v: best(f'{kd}/{v}_{q}_r1/epa_result.jplace') for v in ('stock', 'bug1only', 'bug2only', 'fix', 'stock_rsoff')
             if os.path.exists(f'{kd}/{v}_{q}_r1/epa_result.jplace')}
        s = []
        for a, b in itertools.combinations(B, 2):
            same = sum(B[a][k] == B[b].get(k) for k in B[a]) / len(B[a])
            s.append(f'{a}={b}:{100*same:.0f}%')
        print(os.path.basename(kd), q, ' '.join(s))
