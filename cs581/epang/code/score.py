"""Score jplace files by delta error on the backbone tree.

Usage: python score.py <datadir> <out.tsv> <label=jplace> [<label=jplace> ...]
datadir holds rx.raxml.bestTree (backbone T), backbone.fa, true tree path in true_tree.txt.
Placements on subtrees are mapped back onto T (see deltalib.Backbone.map_jplace).
Appends rows: label, query, delta, lwr
"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deltalib


def load_backbone(d):
    B = [l[1:].strip() for l in open(f'{d}/backbone.fa') if l.startswith('>')]
    true = open(f'{d}/true_tree.txt').read().strip()
    return deltalib.Backbone(f'{d}/rx.raxml.bestTree', true, B)


def main():
    d, out = sys.argv[1], sys.argv[2]
    bb = load_backbone(d)
    with open(out, 'a') as f:
        for arg in sys.argv[3:]:
            lab, path = arg.split('=', 1)
            pl = {}
            for p in path.split(','):
                pl.update(bb.map_jplace(p))
            for q, (c, lwr) in sorted(pl.items()):
                f.write(f'{lab}\t{q}\t{bb.delta(q, c)}\t{lwr:.4f}\n')
            ds = [bb.delta(q, c) for q, (c, _) in pl.items()]
            print(f'{lab}\tn={len(ds)}\tmean_delta={sum(ds)/max(1,len(ds)):.3f}', flush=True)


if __name__ == '__main__':
    main()
