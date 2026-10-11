"""Score all finished nested-experiment runs. Output TSV: center k mode qtype query delta lwr
Usage: python score_nested.py <datadir> <out.tsv> [extra_root ...]"""
import sys, os, json, glob, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score import load_backbone


def main():
    d, out = sys.argv[1], sys.argv[2]
    roots = sys.argv[3:] or ['nested']
    bb = load_backbone(d)
    rows = []
    for root in roots:
        meta = json.load(open(f'{d}/{root}/meta.json'))
        qc = {q: c for c, v in meta['centers'].items() for q in v['queries']}
        for jp in sorted(glob.glob(f'{d}/{root}/*/*/epa_result.jplace') +
                         glob.glob(f'{d}/{root}/*/*/*/epa_result.jplace') +
                         glob.glob(f'{d}/{root}/*/*/*/pplacer.jplace') +
                         glob.glob(f'{d}/{root}/*/*/apples.jplace') +
                         glob.glob(f'{d}/{root}/*/*/*/apples.jplace')):
            parts = jp.split('/')
            mq = parts[-2]
            m = re.match(r'(.+)_(frag|full)$', mq)
            if not m:
                continue
            mode, qt = m.groups()
            if parts[-3] == 'all':
                cen, k = None, 9000
            else:
                cen, k = parts[-4], int(parts[-3][1:])
            try:
                pl = bb.map_jplace(jp)
            except Exception as e:
                print('FAIL', jp, e); continue
            for q, (c, lwr) in pl.items():
                rows.append((cen or qc.get(q, '?'), k, mode, qt, q, bb.delta(q, c), lwr))
    with open(out, 'w') as f:
        f.write('center\tk\tmode\tqtype\tquery\tdelta\tlwr\n')
        for r in rows:
            f.write('\t'.join(map(str, r[:-1])) + f'\t{r[-1]:.4f}\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main()
