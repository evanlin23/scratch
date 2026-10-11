AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R25, SIMHIGH_R26: vote-model tree replicates (helper p25)

Per replicate (AliSim, cs581/gcmtrees/code/gen.sh seeds): one MAGUS draw (paper flags) via
`cs581/gcmtrees/code/run.sh NAME`, then `gg.py run REP 'linsi#es4'`, then vote.py variants via
`cs581/gcmvote/code/run.py` in a separate rep dir /opt/work/gcmvote/reps/NAME (`inputs`, `true.fasta`
symlinked from the gcmtrees rep). The whole lane is `scripts/pipe.sh NAME`.

- `magus.jsonl`: MAGUS end-to-end rows (bbtool_bench).
- `aln.jsonl`: per replicate, gg.py rows (`variant` linsi = MAGUS's merge, recipe = wsoft0.03:linsi&fftns2#es4,
  linsi#es3, linsi&fftns2-op3, linsi#es4; no `B` key) followed by vote.py rows (`B` = 10; magus, hard,
  hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree), 4 FastTree
  processes at a time. Methods `true`, `magus` (gg linsi), `es4` (gg linsi#es4), `recipe`, `vote_hard`, `vote_hard-bb`.
- Bank: `cs581/gcmvote/bank/NAME.tar.gz` + MANIFEST.tsv row (`scripts/bank.sh`), same layout as helper t20's
  bank plus `subsets.json` (from vote/magus).

Checks (`scripts/norm.py`: md5 of row-sorted alignment):
- SIMHIGH_R25: vote `magus` equals gg `linsi` (MAGUS's merge): same row-sorted md5, identical SPFN/SPFP/TC.
  vote `hard` and `hard-bb` keep the same 99881 edges and give the same alignment up to row order (same md5,
  same SP scores); both trees were still built separately (same nRF).
- SIMHIGH_R26: vote `magus` equals gg `linsi` (same row-sorted md5, identical SP scores). vote `hard` (102329 edges)
  and `hard-bb` (100015 edges) differ.
