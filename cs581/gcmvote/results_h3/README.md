AI-assisted (Claude), exploration code for CS581 project

# gcmvote helper h3: fresh MAGUS draws on held-out BAliBASE (BBA0154, BBA0190)

Per dataset, one fresh MAGUS draw with the paper's flags (job line `NAME cs581/data/balibase_clean/RV100_NAME.fasta 25`,
`gcmx.bbtool_bench --draws 0`, then `protcons/code/pc.py rep`), as `gcmgen/code/fresh.sh` does.

- `NAME/magus_draw.jsonl`: the bbtool_bench row for the draw (end-to-end MAGUS score, backbone stats).
- `NAME/baselines.jsonl`: `gcmgen/code/gg.py run REP linsi 'linsi#es3' 'linsi#es4' 'linsi#es5'` (`linsi` = MAGUS's merge).
- `NAME/vote.jsonl`: `gcmvote/code/run.py` (commit 832d0e2 on claude/cs581-gcmvote) on the same rep, every
  PREREG variant at B = 10: magus, es4, hard, hard-bb, soft, soft-bb, soft2, soft4, frac0.2-0.5
  (BBA0154 also has hard+mask; its alignment score equals hard). Scored with gcmx.bbtool_bench.acc_ref
  (FastSP vs the reference), as gg.py does. SP error = avgErr = (SPFN+SPFP)/2.
- `NAME/model_VARIANT.json`: vote.py's fitted mixtures, exposure histogram (n_hist), effective cutoff per n.
- `NAME/edges_hard.npz`: per-edge support k, exposure n, weight, posterior (binomial) from the `hard` run.
