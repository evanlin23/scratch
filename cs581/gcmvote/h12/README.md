AI-assisted (Claude), exploration code for CS581 project

Helper h12: 16S.T (CRW, 5,548 sequences, /opt/data/Datasets/Gutell/16S.T/R0/true_align_clean.txt from the
MAGUS paper's Datasets.zip, Illinois Data Bank doi:10.13012/B2IDB-2643961_V1, fetched by cs581/code/setup.sh).

- magus.sh     one fresh MAGUS draw, paper flags, K = 100 (as the paper ran 16S.T), vectorized graph builder
               (fastbench.py = gcmx.bbtool_bench with --gcmx-fastgraph true), then pc.py rep
- pipeline.sh  gg.py baselines (linsi, linsi#es3/4/5; MAGUS's original graph builder) then every
               pre-declared vote variant from claude/cs581-gcmvote (run.py, B = 10)
- memwatch.py  wall time, peak summed RSS of the process tree, 4 h limit per step
Rows: ../results_h12/aln.jsonl (in progress).

## What happened (2026-10-11, 4 cores, 15 GB)

- MAGUS draw: 4,386 s wall (15,189 s CPU), peak RSS 3.5 GB (whole bench incl. its merge-mafft control, 494 s).
- Graph check: gg.py `linsi` (MAGUS's original builder) gives the same alignment as the draw up to column order
  (identical column multiset and FastSP); vote.py `magus` also reproduces it.
- Merges: `linsi` alone; the other 15 three at a time on the shared 4 cores (`concurrent_jobs` in each row), so
  those walls are not single-job walls. Peak RSS ≤ 3.9 GB per merge. Nothing hit the 4 h limit.
- `hard+mask` not run (its scored alignment equals `hard`; skipped per the orchestrator). B = 5/20 not run.
- Rows: ../results_h12/aln.jsonl (err_pct = (SPFN+SPFP)/2 × 100). Rep bank: ../bank/16S.T_R0.tar.gz.
