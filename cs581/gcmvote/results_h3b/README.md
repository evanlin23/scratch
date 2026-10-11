AI-assisted (Claude), exploration code for CS581 project

# gcmvote helper h3b: BBA0081, BBA0117 (held-out)

One fresh MAGUS draw per set (paper flags, as gcmgen/code/fresh.sh; `fresh_magus.jsonl`), rep dir via `pc.py rep`,
reference = cs581/data/balibase_clean/RV100_NAME.fasta. Run: `bash pipeline.sh NAME`.

- `NAME.gg_results.jsonl`: `gg.py run REP linsi linsi#es3 linsi#es4 linsi#es5`
- `NAME.results.jsonl`: vote.py variants (B = 10), scored as run.py/gg.py (`acc_ref` = FastSP vs the reference)
- `NAME.VARIANT.model.json`: vote.py's model.json (fit, cutoff_by_n, support/exposure histograms)

SP error % = avgErr × 100:

| set | MAGUS (linsi) | es3 | es4 | es5 | vote magus | vote es4 | hard | soft | soft2 | soft4 | hard-bb | soft-bb |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BBA0081 | 58.84 | 58.71 | 57.85 | 56.97 | 58.84 | 58.05 | 57.13 | 56.07 | 55.83 | 55.18 | 63.06 | 56.46 |
| BBA0117 | 13.35 | 13.31 | 13.31 | 13.31 | 13.35 | 13.31 | 12.92 | 12.92 | 12.92 | 12.92 | 12.91 | 12.93 |

vote.py `magus` reproduces gg.py `linsi` exactly on both sets.
