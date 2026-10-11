AI-assisted (Claude), exploration code for CS581 project

# Helper h5d: RNASim_R1 (held-out)

One fresh MAGUS draw (paper flags; reference /opt/data/Datasets/RNASim/1000/R1/true_align.txt), built and run with `pipeline.sh`.

- `RNASim_R1.fresh_magus.jsonl`: bbtool_bench rows for the MAGUS draw
- `RNASim_R1.gg_results.jsonl`: gg.py rows (linsi, linsi#es3, linsi#es4, linsi#es5)
- `RNASim_R1.results.jsonl`: vote.py rows scored by gcmvote/code/run.py (FastSP vs the reference; same fields as gg.py plus
  B, kept_edges, edges, kept_weight_frac, cutoff_by_n)
- `RNASim_R1.<variant>.model.json`: vote.py mixture fits, exposure histogram (n_hist), cutoffs
- Check: vote.py `magus` out.fasta has the same sequences, row for row, as gg.py `linsi` (only the sequence order differs); the scores are identical.

| variant | SP error % | SPFN % | SPFP % | TC |
|---|---|---|---|---|
| MAGUS (gg linsi) | 8.95 | 9.10 | 8.80 | 0.0459 |
| linsi#es3 | 8.93 | 9.01 | 8.84 | 0.0471 |
| linsi#es4 | 8.88 | 8.96 | 8.80 | 0.0473 |
| linsi#es5 | 8.92 | 9.02 | 8.82 | 0.0473 |
| vote magus | 8.95 | 9.10 | 8.80 | 0.0459 |
| vote es4 | 8.88 | 8.96 | 8.80 | 0.0473 |
| hard | 9.23 | 10.38 | 8.08 | 0.0471 |
| soft | 9.13 | 10.13 | 8.13 | 0.0453 |
| soft2 | 9.22 | 10.39 | 8.06 | 0.0459 |
| soft4 | 9.25 | 10.46 | 8.04 | 0.0459 |
| hard-bb | 9.22 | 10.37 | 8.08 | 0.0473 |
| soft-bb | 9.16 | 10.20 | 8.12 | 0.0447 |
