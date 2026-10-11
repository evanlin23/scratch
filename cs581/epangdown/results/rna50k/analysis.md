queries placed by all runs: 1000 of 1000

| run | mean delta | median | % delta=0 | vs stock_b2000_frag_s0: diff | W/T/L | Wilcoxon p | wall (m:s) | max RSS (GB) |
|---|---|---|---|---|---|---|---|---|
| fix_b2000_frag_s0 | 0.530 | 0 | 68 | +0.000 | 0/1000/0 | 1 | 2:57.34 | 2.8 |
| stock_b2000_frag_s0 | 0.530 | 0 | 68 |  |  |  | 3:05.73 | 2.8 |
| fix_b5000_frag_s0 | 0.507 | 0 | 68 | -0.023 | 22/963/15 | 0.089 | 5:11.32 | 7.4 |
| stock_b5000_frag_s0 | 1.434 | 1 | 41 | +0.904 | 78/488/434 | 3.9e-53 | 6:30.41 | 7.4 |

fix vs stock at the same subtree size (W = fix lower):
  fix_b2000_frag_s0 vs stock_b2000_frag_s0: diff/WTL/p = ('+0.000', '0/1000/0', '1')
  fix_b5000_frag_s0 vs stock_b5000_frag_s0: diff/WTL/p = ('-0.927', '439/487/74', '1e-55')
