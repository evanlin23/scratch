queries placed by all runs: 1000 of 1000

| run | mean delta | median | % delta=0 | vs stock_b2000_frag: diff | W/T/L | Wilcoxon p | wall (m:s) | max RSS (GB) |
|---|---|---|---|---|---|---|---|---|
| stock_b2000_frag | 0.811 | 0 | 62 |  |  |  | 0:43.17 | 2.9 |
| fix_whole_frag | 0.777 | 0 | 63 | -0.034 | 40/938/22 | 0.068 |  | nan |
| stock_whole_frag | 1.719 | 1 | 35 | +0.908 | 103/430/467 | 1.9e-48 |  | nan |

fix vs stock at the same subtree size (W = fix lower):
  fix_whole_frag vs stock_whole_frag: diff/WTL/p = ('-0.942', '470/436/94', '5.7e-52')
