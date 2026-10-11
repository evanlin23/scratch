AI-assisted (Claude), exploration code for CS581 project

# gcmvote helper h11: 16S.3 (CRW, 5,489 seqs, clean reference)

Data: /opt/data/Datasets/Gutell/16S.3/R0/true_align_clean.txt from the MAGUS paper's Datasets.zip
(Illinois Data Bank doi:10.13012/B2IDB-2643961_V1, via cs581/code/setup.sh). Machine: 4 cores, 15 GB RAM.

- MAGUS draw: `GCMX_FASTGRAPH=true python3 -m gcmx.bbtool_bench` (job `16S.3_R0 ... 100`, paper flags,
  K = 100 as in the paper's 16S.3 run), then `pc.py rep`. The vectorized builder is the documented
  identical-graph drop-in; the merge-only control with it reproduced MAGUS exactly (avgErr 10.170).
  Wall 4,763 s (1:28:44), CPU 16,430 s, peak RSS 3.40 GB (/usr/bin/time -v, largest process).
- Baselines: `gg.py run REP linsi 'linsi#es3' 'linsi#es4' 'linsi#es5'` (harness unchanged: pure-Python graph).
- Vote: `run.py REP magus es4 hard hard-bb soft soft-bb soft2 soft4 frac0.2 frac0.3 frac0.4 frac0.5` from
  branch claude/cs581-gcmvote @832d0e2 (B = 10). hard+mask not scored (orchestrator: alignment equals hard).
- Merges ran 2-3 at a time with -np 4 each, so merge_wall values are under contention. Peak RSS per lane
  (largest process): gg lanes 4.05 / 3.61 GB, vote lanes 3.43 / 3.43 GB.
- aln.jsonl: one row per method; avgErr = (SPFN + SPFP) / 2 from FastSP (gcmx.bbtool_bench.acc_ref).
