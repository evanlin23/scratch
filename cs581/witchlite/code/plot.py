#!/usr/bin/env python3
"""Accuracy-vs-runtime plot: dSPFN vs WITCH (pp) against % of WITCH wall time, one panel per instance.
usage: plot.py RESULTS_DIR OUT.png"""
import json, sys, glob, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GROUPS = [  # (label, color, marker, method keys)
    ('WITCH top-k / adaptive k (all HMMs scored)', '#2a78d6', 'o', ['all/k1', 'all/k3', 'all/k10_t99', 'all/k10_t95']),
    ('hier descent / EarlyStop (UPP2-style)', '#eb6834', 's', ['hier/k10', 'hier/k10_t99', 'hier_es/k10']),
    ('beam-b descent (b=2,3,4)', '#1baf7a', 'D', ['beam/k10', 'beam/k10_t99', 'beam3/k10', 'beam3/k10_t99', 'beam4/k10']),
    ('hybrid BLAST/beam-3', '#e87ba4', 'P', ['hyb/b3_100']),
    ('BLAST-guided path (+siblings)', '#6250d6', '^', ['blastpath/k10', 'blastpath/k10_t99', 'blastpath_sib/k10',
                                                        'blastpath_sib/k10_t99', 'blastpath_sib/k10_t95']),
    ('BLASTN only (-task blastn)', '#52514e', 'x', ['BLAST-sens']),
]
files = [f for f in sorted(glob.glob(f'{sys.argv[1]}/*.json')) if 'methods' in json.load(open(f))]
fig, axes = plt.subplots(1, len(files), figsize=(4.2 * len(files), 3.8), squeeze=False)
for ax, f in zip(axes[0], files):
    r = json.load(open(f)); M = r['methods']; wall = r['witch_default_time']['wall']
    base = M['WITCH']['spfn']
    for lab, col, mk, keys in GROUPS:
        xs = [100 * M[k]['time'] / wall for k in keys if k in M]
        ys = [100 * (M[k]['spfn'] - base) for k in keys if k in M]
        ax.scatter(xs, ys, c=col, marker=mk, s=40, label=lab, zorder=3)
    ax.scatter([100], [0], c='#0b0b0b', marker='*', s=90, label='WITCH default', zorder=3)
    ax.axhline(1.0, color='#999', lw=1, ls='--'); ax.axvline(20, color='#999', lw=1, ls='--')
    ax.set_title(f"{os.path.basename(f)[:-5]} (WITCH SPFN {100 * base:.2f}%)", fontsize=10)
    ax.set_xlabel('% of WITCH wall time'); ax.set_ylabel('SPFN − WITCH (pp)')
    ax.set_xlim(-3, 105); top = max(2.0, min(10, max(100 * (m['spfn'] - base) for k, m in M.items() if k != 'BLAST') + 0.3))
    bot = min(-0.3, min(100 * (m['spfn'] - base) for m in M.values()) - 0.2)
    ax.set_ylim(bot, top); ax.grid(alpha=0.25)
axes[0][0].legend(fontsize=7, loc='upper right')
fig.suptitle('Query-to-backbone SPFN vs runtime (go/kill box: < 20% time, < 1 pp SPFN)', fontsize=10)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=130)
