"""Score several jplace files by delta error on the backbone (memory-lean version of
cs581/epang/code/score.py: per-query caches are dropped after each query, needed for 77K-leaf
backbones). Usage: python score_multi.py <datadir> <out.tsv> <label=jplace> [...]
Writes rows: label, query, delta, lwr (overwrites out.tsv)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deltalib
from score import load_backbone


def main():
    d, out = sys.argv[1], sys.argv[2]
    bb = load_backbone(d)
    pls = []
    for arg in sys.argv[3:]:
        lab, path = arg.split('=', 1)
        pls.append((lab, bb.map_jplace(path)))
    qs = sorted(set().union(*[set(p) for _, p in pls]))
    with open(out, 'w') as f:
        f.write('label\tquery\tdelta\tlwr\n')
        for q in qs:
            for lab, p in pls:
                if q in p:
                    c, lwr = p[q]
                    f.write(f'{lab}\t{q}\t{bb.delta(q, c)}\t{lwr:.4f}\n')
            bb._q.clear()
    print('scored', len(qs), 'queries x', len(pls), 'runs', flush=True)


if __name__ == '__main__':
    main()
