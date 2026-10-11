queries placed by all runs: 1000 of 1000

| run | mean delta | median | % delta=0 | vs stock_b2000_frag: diff | W/T/L | Wilcoxon p | wall (m:s) | max RSS (GB) |
|---|---|---|---|---|---|---|---|---|
| stock_b2000_frag | 0.791 | 0 | 58 |  |  |  | 1:08.80 | 2.5 |
| fix_b5000_frag | 0.796 | 0 | 58 | +0.005 | 13/970/17 | 0.48 | 1:22.79 | 6.6 |
| fix_whole_frag | 0.776 | 0 | 59 | -0.015 | 14/977/9 | 0.25 |  | nan |
| stock_whole_frag | 2.549 | 2 | 25 | +1.758 | 60/379/561 | 1.7e-80 |  | nan |

fix vs stock at the same subtree size (W = fix lower):
  fix_whole_frag vs stock_whole_frag: diff/WTL/p = ('-1.773', '565/381/54', '2.1e-82')
