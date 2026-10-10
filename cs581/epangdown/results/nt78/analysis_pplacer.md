queries placed by all runs: 1000 of 1000

| run | mean delta | median | % delta=0 | vs pplacer_b2000_frag: diff | W/T/L | Wilcoxon p | wall (m:s) | max RSS (GB) |
|---|---|---|---|---|---|---|---|---|
| fix_b2000_frag | 1.721 | 1 | 37 | +0.016 | 77/818/105 | 0.27 | 5:36.30 | 2.4 |
| pplacer_b2000_frag | 1.705 | 1 | 41 |  |  |  | 13:30.04 | 2.0 |
| stock_b2000_frag | 1.720 | 1 | 37 | +0.015 | 73/823/104 | 0.16 | 5:05.86 | 2.4 |
| fix_b5000_frag | 1.697 | 1 | 37 | -0.008 | 89/789/122 | 0.33 | 7:35.79 | 6.4 |
| pplacer_b5000_frag | 1.633 | 1 | 42 | -0.072 | 46/923/31 | 0.041 | 23:31.12 | 5.1 |
| stock_b5000_frag | 4.769 | 3 | 16 | +3.064 | 107/287/606 | 3e-80 | 9:15.77 | 6.5 |
| fix_b10000_frag | 1.690 | 1 | 38 | -0.015 | 99/776/125 | 0.72 | 8:15.05 | 13.0 |

fix vs stock at the same subtree size (W = fix lower):
  fix_b2000_frag vs stock_b2000_frag: diff/WTL/p = ('+0.001', '10/981/9', '0.81')
  fix_b5000_frag vs stock_b5000_frag: diff/WTL/p = ('-3.072', '584/325/91', '9.7e-81')
