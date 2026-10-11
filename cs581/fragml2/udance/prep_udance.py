#!/usr/bin/env python3
"""Prepare a uDance single-gene workdir: alignments/gene.fasta (all 1000 seqs, true alignment),
backbone.nwk (FastTree tree on full-length seqs, ungapped len >= 0.5*median)."""
import sys, os, re, statistics, subprocess
rep_dir, wd = sys.argv[1], sys.argv[2]
aln = os.path.join(rep_dir, 'true_align.fasta')
seqs = {}; name = None
for line in open(aln):
    line = line.strip()
    if line.startswith('>'): name = line[1:].split()[0]; seqs[name] = []
    elif name: seqs[name].append(line)
seqs = {k: ''.join(v).upper() for k, v in seqs.items()}
ul = {k: len(v.replace('-', '')) for k, v in seqs.items()}
med = statistics.median(ul.values())
bb = sorted(k for k in seqs if ul[k] >= 0.5 * med)
os.makedirs(os.path.join(wd, 'alignments'), exist_ok=True)
with open(os.path.join(wd, 'alignments', 'gene.fasta'), 'w') as f:
    for k, v in seqs.items(): f.write('>%s\n%s\n' % (k, v))
tre = os.path.join(rep_dir, 'trees/cache/true_align.bb_ft_0.5.tre')
if os.path.exists(tre):
    t = open(tre).read().strip()
else:
    bbfa = os.path.join(wd, 'bb.fasta')
    with open(bbfa, 'w') as f:
        for k in bb: f.write('>%s\n%s\n' % (k, seqs[k]))
    t = subprocess.run('nice -n 10 /usr/bin/FastTree -nt -gtr -gamma -quiet < %s' % bbfa, shell=True,
                       capture_output=True, text=True, check=True).stdout.strip()
labels = set(re.findall(r'[(,]([^(),:;]+)', t))
print('median ungapped', med, 'backbone', len(bb), 'tree leaves', len(labels), 'match', labels == set(bb),
      'min query len', min(ul[k] for k in seqs if k not in labels), file=sys.stderr)
open(os.path.join(wd, 'backbone.nwk'), 'w').write(t + '\n')
