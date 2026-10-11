AI-assisted (Claude), exploration code for CS581 project

Helper h5b: held-out ROSE 1000M2 R1 and 1000L1 R1 (one fresh MAGUS draw each, paper flags; run with
cs581/gcmvote/code_h5b/pipe.sh, which follows gcmgen/code/fresh.sh).
- NAME.gg.jsonl: gg.py baselines (linsi = MAGUS, linsi#es3/4/5)
- NAME.results.jsonl: vote.py variants (magus es4 hard soft soft2 soft4 hard-bb soft-bb), B = 10, run.py's row format,
  scored with code_h5b/par.py (run.py's scoring, variants run in parallel)
- NAME.VARIANT.model.json: vote.py's model.json; NAME.calib.json: calib.py on the magus edges
- fresh_magus.jsonl: bbtool_bench output of the two draws
vote.py `magus` reproduces gg.py `linsi` exactly on both replicates (same SPFN, SPFP, TC).
