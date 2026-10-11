# Frozen choice (from TRAIN only, by the PREREG.md rule), written before the held-out run of the chosen variant

Mean FN rate on TRAIN (n = 40: 36 FastMulRFS + 4 DISCO):

| id | flags | mean FN |
|---|---|---|
| U | none | 0.07162 |
| **S** | `-W sup` | **0.07059** |
| SE | `-W supe` | 0.0719 |
| S2 | `-W sup -p 2` | 0.0716 |
| L | `-W len -N diam` | 0.0801 |
| C1 / C3 / C5 | `-c 0.1` / `-c 0.3` / `-c 0.5` | 0.0745 / 0.0745 / 0.0721 |
| N | `-N mean` | 0.0734 |
| P | `-A pair` | 0.0783 |
| G | `-G msup` | 0.07136 (the only lever that beat U on its own) |
| S+G | `-W sup -G msup` | 0.07136 (worse than S, so G is not kept) |

**Steps of the rule.**
1. Best weighting: S.
2. Levers that helped alone: G only. Adding it to S raises FN, so it is dropped.
3. S beats U by 0.00103, which is outside the 0.001 tie band.

**Chosen: wASTRID-Pro = `wapro -M pro -u -W sup` + FastME 2.**

**Baseline setting chosen on TRAIN.**
- Asteroid: `--use-gene-bl` (0.0721) beats the default (0.0732), so `asteroid-bl` is the Asteroid baseline on
  HELD-OUT.
- Other baselines use their pre-registered settings.

**Compute note (pre-registered drop rule applied).**
- On the 1000-gene DISCO data, wQFM-GDL takes about 35–45 min, DISCO+wASTRAL about 25 min and ASTRAL-Pro3 about
  13 min per dataset on this 4-core machine.
- The fast methods run on all 12 DISCO held-out datasets.
- The three slow baselines run on DISCO held-out datasets in this priority order: `default` 03, 04;
  `gdl_1e-10_05` 03, 04; `gdl_5e-10_05` 01–04; `ils_2e8` 01–04, until time runs out.
- Paired tests use the intersection of datasets; coverage is reported.
