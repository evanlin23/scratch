queries placed by all runs: 1000 of 1000

| run | mean delta | median | % delta=0 | vs stock_b2000_frag: diff | W/T/L | Wilcoxon p | wall (m:s) | max RSS (GB) |
|---|---|---|---|---|---|---|---|---|
| fix_b2000_frag | 25.439 | 13 | 8 | +0.167 | 306/414/280 | 0.56 | 2:01.24 | 2.8 |
| stock_b2000_frag | 25.272 | 14 | 8 |  |  |  | 2:44.49 | 2.8 |
| bug1only_b5000_frag | 28.203 | 17 | 4 | +2.931 | 391/115/494 | 0.00027 | 2:24.46 | 7.4 |
| bug2only_b5000_frag | 25.330 | 13 | 8 | +0.058 | 335/351/314 | 0.58 | 3:14.96 | 7.4 |
| fix_b5000_frag | 26.850 | 13 | 8 | +1.578 | 323/355/322 | 0.45 | 2:32.55 | 7.4 |
| stock_b5000_frag | 28.026 | 17 | 4 | +2.754 | 389/125/486 | 0.00055 | 4:11.60 | 7.4 |

fix vs stock at the same subtree size (W = fix lower):
  fix_b2000_frag vs stock_b2000_frag: diff/WTL/p = ('+0.167', '306/414/280', '0.56')
  fix_b5000_frag vs stock_b5000_frag: diff/WTL/p = ('-1.176', '485/122/393', '0.0092')
  bug2only_b5000_frag vs fix_b5000_frag: diff/WTL/p = ('-1.520', '323/389/288', '0.11'); identical delta on 38.9% of queries
  bug2only_b5000_frag vs stock_b5000_frag: diff/WTL/p = ('-2.696', '482/131/387', '0.00065'); identical delta on 13.1% of queries
  bug1only_b5000_frag vs fix_b5000_frag: diff/WTL/p = ('+1.353', '393/130/477', '0.0059'); identical delta on 13.0% of queries
  bug1only_b5000_frag vs stock_b5000_frag: diff/WTL/p = ('+0.177', '69/874/57', '0.71'); identical delta on 87.4% of queries
