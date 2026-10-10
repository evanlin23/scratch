| regime | k | S structure | claw | N | true_found | continuum | alt_k'<k | alt_k'=k | any_alt_k'<=k | alt_k'=k+1 |
|---|---|---|---|---|---|---|---|---|---|---|
| generic | 1 | matching | noclaw | 1279 | 1279 | 0 | 0 | 0 | 0 | 0 |
| generic | 2 | adjacent | noclaw | 1498 | 1498 | 1498 | 0 | 0 | 0 | 1498 |
| generic | 2 | matching | noclaw | 2997 | 2997 | 0 | 0 | 0 | 0 | 1506 |
| generic | 3 | adjacent | claw | 1043 | 1043 | 1043 | 782 | 743 | 799 | 1043 |
| generic | 3 | adjacent | noclaw | 2893 | 2893 | 2893 | 787 | 2006 | 2006 | 2893 |
| generic | 3 | matching | claw | 591 | 591 | 530 | 0 | 167 | 167 | 591 |
| generic | 3 | matching | noclaw | 598 | 598 | 0 | 0 | 0 | 0 | 583 |
| generic_y0 | 1 | matching | noclaw | 1279 | 1279 | 0 | 0 | 0 | 0 | 0 |
| generic_y0 | 2 | adjacent | noclaw | 1498 | 1498 | 1498 | 0 | 0 | 0 | 0 |
| generic_y0 | 2 | matching | noclaw | 2997 | 2997 | 0 | 0 | 0 | 0 | 1506 |
| generic_y0 | 3 | adjacent | claw | 1043 | 1043 | 1043 | 789 | 755 | 804 | 1043 |
| generic_y0 | 3 | adjacent | noclaw | 2893 | 2893 | 2893 | 787 | 1172 | 1959 | 1291 |
| generic_y0 | 3 | matching | claw | 591 | 591 | 530 | 0 | 164 | 164 | 591 |
| generic_y0 | 3 | matching | noclaw | 598 | 598 | 0 | 0 | 0 | 0 | 583 |
| sym | 1 | matching | noclaw | 1279 | 1279 | 0 | 0 | 0 | 0 | 0 |
| sym | 2 | adjacent | noclaw | 1498 | 1498 | 1498 | 0 | 0 | 0 | 1498 |
| sym | 2 | matching | noclaw | 2997 | 2997 | 0 | 0 | 0 | 0 | 1506 |
| sym | 3 | adjacent | claw | 1043 | 1043 | 1043 | 0 | 260 | 260 | 1043 |
| sym | 3 | adjacent | noclaw | 2893 | 2893 | 2893 | 787 | 1808 | 1808 | 2893 |
| sym | 3 | matching | claw | 591 | 591 | 530 | 0 | 0 | 0 | 591 |
| sym | 3 | matching | noclaw | 598 | 598 | 0 | 0 | 0 | 0 | 583 |
| sym_y0 | 1 | matching | noclaw | 1279 | 1279 | 0 | 0 | 0 | 0 | 0 |
| sym_y0 | 2 | adjacent | noclaw | 1498 | 1498 | 1498 | 0 | 0 | 0 | 0 |
| sym_y0 | 2 | matching | noclaw | 2997 | 2997 | 0 | 0 | 0 | 0 | 1506 |
| sym_y0 | 3 | adjacent | claw | 1043 | 1043 | 1043 | 0 | 260 | 260 | 1043 |
| sym_y0 | 3 | adjacent | noclaw | 2893 | 2893 | 2893 | 787 | 1021 | 1808 | 1527 |
| sym_y0 | 3 | matching | claw | 591 | 591 | 530 | 0 | 0 | 0 | 591 |
| sym_y0 | 3 | matching | noclaw | 598 | 598 | 0 | 0 | 0 | 0 | 583 |
| ultra | 1 | matching | noclaw | 1279 | 1279 | 0 | 0 | 0 | 0 | 0 |
| ultra | 2 | adjacent | noclaw | 1498 | 1498 | 1498 | 0 | 0 | 0 | 973 |
| ultra | 2 | matching | noclaw | 2997 | 2997 | 0 | 0 | 0 | 0 | 1506 |
| ultra | 3 | adjacent | claw | 1043 | 1043 | 1043 | 684 | 638 | 701 | 1043 |
| ultra | 3 | adjacent | noclaw | 2893 | 2893 | 2893 | 787 | 1600 | 1876 | 2410 |
| ultra | 3 | matching | claw | 591 | 591 | 530 | 0 | 134 | 134 | 591 |
| ultra | 3 | matching | noclaw | 598 | 598 | 0 | 0 | 0 | 0 | 583 |

k=3, S a matching: (regime, claw, has alternative with k'<=k) -> count
  ('generic', False, False): 598
  ('generic', True, False): 424
  ('generic', True, True): 167
  ('generic_y0', False, False): 598
  ('generic_y0', True, False): 427
  ('generic_y0', True, True): 164
  ('sym', False, False): 598
  ('sym', True, False): 591
  ('sym_y0', False, False): 598
  ('sym_y0', True, False): 591
  ('ultra', False, False): 598
  ('ultra', True, False): 457
  ('ultra', True, True): 134

examples (generic regime, S a matching, alternative with k'<=k):
  {"n": 5, "newick": "(2:1,(0:1,3:1)i6:1,(1:1,4:1)i7:1)i5;", "S": [[2, 5], [0, 6], [1, 7]], "alts": [[[0, 6], [6, 5], [1, 7]], [[0, 6], [1, 7], [7, 5]]]}
  {"n": 5, "newick": "(2:1,(0:1,3:1)i6:1,(1:1,4:1)i7:1)i5;", "S": [[2, 5], [0, 6], [7, 4]], "alts": [[[0, 6], [6, 5], [7, 4]], [[0, 6], [7, 5], [7, 4]]]}
  {"n": 5, "newick": "(2:1,(0:1,3:1)i6:1,(1:1,4:1)i7:1)i5;", "S": [[2, 5], [6, 3], [1, 7]], "alts": [[[6, 5], [6, 3], [1, 7]], [[6, 3], [1, 7], [7, 5]]]}
  {"n": 5, "newick": "(2:1,(0:1,3:1)i6:1,(1:1,4:1)i7:1)i5;", "S": [[2, 5], [6, 3], [7, 4]], "alts": [[[6, 5], [6, 3], [7, 4]], [[6, 3], [7, 5], [7, 4]]]}
  {"n": 5, "newick": "(1:1,(0:1,3:1)i6:1,(2:1,4:1)i7:1)i5;", "S": [[1, 5], [0, 6], [2, 7]], "alts": [[[0, 6], [6, 5], [2, 7]], [[0, 6], [2, 7], [7, 5]]]}
  {"n": 5, "newick": "(1:1,(0:1,3:1)i6:1,(2:1,4:1)i7:1)i5;", "S": [[1, 5], [6, 3], [7, 4]], "alts": [[[6, 5], [6, 3], [7, 4]], [[6, 3], [7, 5], [7, 4]]]}
