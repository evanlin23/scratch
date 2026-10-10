#!/bin/bash
# lnL (GTR+G, raxml-ng --evaluate) of the true tree and Park et al.'s published trees on RNASim10K R2:
# does the time-capped RAxML-NG search leave likelihood on the table at 10K taxa (search headroom)?
cd "$(dirname "$0")/.."
D=/opt/data/mlcache/ParkRNASim10K/R2
A=$D/true_align.fasta.clean.fasta
for t in true_tree published_gtm published_raxml published_iq published_ft; do
  python3 -c "import sys,json; sys.path.insert(0,'code'); import evaltree, treeerr; print(json.dumps({'tree':'$t','lnl':evaltree.evaluate('$A','$D/$t.tre'),'fn':treeerr.error('$D/true_tree.tre','$D/$t.tre')['fn_rate']}))" >> results/eval10k.jsonl &
done
wait
