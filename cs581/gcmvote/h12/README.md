AI-assisted (Claude), exploration code for CS581 project

Helper h12: 16S.T (CRW, 5,548 sequences, /opt/data/Datasets/Gutell/16S.T/R0/true_align_clean.txt from the
MAGUS paper's Datasets.zip, Illinois Data Bank doi:10.13012/B2IDB-2643961_V1, fetched by cs581/code/setup.sh).

- magus.sh     one fresh MAGUS draw, paper flags, K = 100 (as the paper ran 16S.T), vectorized graph builder
               (fastbench.py = gcmx.bbtool_bench with --gcmx-fastgraph true), then pc.py rep
- pipeline.sh  gg.py baselines (linsi, linsi#es3/4/5; MAGUS's original graph builder) then every
               pre-declared vote variant from claude/cs581-gcmvote (run.py, B = 10)
- memwatch.py  wall time, peak summed RSS of the process tree, 4 h limit per step
Rows: ../results_h12/aln.jsonl (in progress).
