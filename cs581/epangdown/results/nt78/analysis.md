queries placed by all runs: 1000 of 1000

| run | mean delta | median | % delta=0 | vs stock_b2000_frag: diff | W/T/L | Wilcoxon p | wall (m:s) | max RSS (GB) |
|---|---|---|---|---|---|---|---|---|
| fix_b2000_frag | 1.721 | 1 | 37 | +0.001 | 10/981/9 | 0.81 | 5:36.30 | 2.4 |
| stock_b2000_frag | 1.720 | 1 | 37 |  |  |  | 5:05.86 | 2.4 |
| bug1only_b5000_frag | 4.756 | 3 | 16 | +3.036 | 103/310/587 | 4.1e-79 | 6:52.66 | 6.5 |
| bug2only_b5000_frag | 1.686 | 1 | 37 | -0.034 | 38/928/34 | 0.37 | 8:59.67 | 6.4 |
| fix_b5000_frag | 1.697 | 1 | 37 | -0.023 | 39/924/37 | 0.56 | 7:35.79 | 6.4 |
| stock_b5000_frag | 4.769 | 3 | 16 | +3.049 | 103/309/588 | 3.1e-79 | 9:15.77 | 6.5 |
| fix_b10000_frag | 1.690 | 1 | 38 | -0.030 | 46/913/41 | 0.32 | 8:15.05 | 13.0 |
| stock_b10000_frag | 5.089 | 3 | 16 | +3.369 | 101/288/611 | 2.7e-84 | 10:38.01 | 13.1 |

fix vs stock at the same subtree size (W = fix lower):
  fix_b2000_frag vs stock_b2000_frag: diff/WTL/p = ('+0.001', '10/981/9', '0.81')
  fix_b5000_frag vs stock_b5000_frag: diff/WTL/p = ('-3.072', '584/325/91', '9.7e-81')
  fix_b10000_frag vs stock_b10000_frag: diff/WTL/p = ('-3.399', '620/286/94', '3e-87')
  bug2only_b5000_frag vs fix_b5000_frag: diff/WTL/p = ('-0.011', '10/983/7', '0.43'); identical delta on 98.3% of queries
  bug2only_b5000_frag vs stock_b5000_frag: diff/WTL/p = ('-3.083', '586/323/91', '3e-81'); identical delta on 32.3% of queries
  bug1only_b5000_frag vs fix_b5000_frag: diff/WTL/p = ('+3.059', '91/326/583', '1.7e-80'); identical delta on 32.6% of queries
  bug1only_b5000_frag vs stock_b5000_frag: diff/WTL/p = ('-0.013', '9/986/5', '0.53'); identical delta on 98.6% of queries
