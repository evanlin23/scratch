# Interim look (deviation from PREREG, requested by the orchestrating session at 23:03 UTC)

**What was run.** `wapro:-W sup` (the best weighting on the FastMulRFS TRAIN data) on the 9 held-out FastMulRFS
datasets that already had all baselines. I sent the result to the orchestrator as an interim number. This happened
*before* `CHOSEN.md` was written, so these 9 datasets (FastMulRFS rep 04, 9 of its 12 condition × length cells) were
seen once.

**What did not change.** The selection rule in PREREG.md is applied to TRAIN only, as written. The final held-out
analysis is reported both with and without these 9 datasets.

wASTRID-Pro (-W sup) minus each method, on n = 9:

| vs | mean diff | W/T/L | p |
|---|---|---|---|
| wQFM-GDL | +0.0080 | 3/2/4 | 0.30 |
| Asteroid | +0.0057 | 3/1/5 | 0.68 |
| ASTRID-DISCO | +0.0034 | 1/4/4 | 0.50 |
| ASTRID-Pro | +0.0023 | 3/2/4 | 0.84 |
| ASTRAL-Pro3 | −0.0115 | 5/2/2 | 0.12 |
