# One MAGUS variant for DNA, RNA and proteins? Self-soft × consensus evidence

Exploration pilot, CS581 (UIUC), Fall 2026, about 7 h on one 4-core machine. Merge-only paired comparisons
against MAGUS on the same subsets and backbone sets. Code: `code/` (`mg.py` driver, `run_all.sh` queue,
`summarize.py` tables). Raw rows: `results/*.results.jsonl`. All tables: `results/tables.md`.

## TL;DR

| variant (no data-type switch, no reference) | DNA/RNA (n = 6) | protein (n = 9) | pooled (n = 15) | extra wall (median) |
|---|---|---|---|---|
| `ss:wsoft0.03:linsi&fftns2`: self-soft on gcmgen's weighted soft consensus | **−0.88, 6/0/0, p = 0.031** | −0.27, 5/1/3, p = 0.46 (1 timeout) | **−0.53, 11/1/3, p = 0.068** | +250 s (20%) |
| `ss:linsi`: self-soft (the published recipe) | −0.74, 6/0/0, p = 0.031 | **+1.28**, 4/1/4, p = 0.95 (1 timeout, one +11.07) | +0.41, 10/1/4, p = 0.15 | +169 s (23%) |
| `wsoft0.03:linsi&fftns2` (gcmgen recipe alone) | −0.08, 2/3/1 | −0.32, 6/0/3 | −0.22, 8/3/4, p = 0.36 | +42 s (4%) |
| `linsi&fftns2-op3`: hard consensus (B) | **+4.50**, 1/1/4 | **−1.31, 8/0/1, p = 0.012** | +1.01, 9/1/5 | ≈ 0 |
| `linsi\|cons0.7`: self-consistency mask | **+3.69**, 2/0/4 | −0.99, 6/2/1, p = 0.027 | +0.88, 8/2/5 | ≈ 0 |
| `soft0.5:linsi&fftns2-op3` | −0.07, 2/4/0 | +0.31, 5/2/2 (BBA0134 +5.08) | +0.16, 7/6/2 | +19 s (3%) |

Δ = error(variant) − error(MAGUS) in points, error = (SPFN + SPFP)/2. W/T/L uses |Δ| < 0.05 as a tie, and
p is a two-sided Wilcoxon test. A self-soft merge that ran longer than 900 s counts as a loss and is left
out of the mean.

**Verdict: not promising as "one variant for everything" (kill under the go/kill rule).**
- Everything with self-soft wins on DNA/RNA (6/0/0).
- Everything with consensus filtering wins on proteins.
- No variant wins significantly on both. The best compromise, self-soft on weighted soft consensus, is
  −0.53 pooled (p = 0.068).
- It loses on 2 of the 3 held-out proteins run after it was chosen (BBA0117 +1.02, SIMMOD_R1 +2.83) and
  times out on SIMHIGH_R1.
- The obstacle is not the evidence. On hard proteins, self-soft's 75-group minclusters trace blows up in
  runtime and in SPFN.

## 1. Question and candidates

Is there one MAGUS variant, with no data-type switch and no reference, that is at least as accurate as MAGUS
on DNA, RNA and proteins and clearly better on average, at modest cost? It combines two earlier pilots:

- **(A) self-soft** (`cs581/code/gcmx/e2e_bench.py`):
  1. Take MAGUS's output restricted to 8 random sequences per subset, giving 10 pseudo-backbones.
  2. HMM-extend them to all sequences (`gcmx.extend`).
  3. Split each subset alignment into 3 similarity groups (`gcmx.split`).
  4. Merge all 75 groups at once with the pseudo-backbones as the only evidence (MAGUS default merge
     flags, MCL with 4 threads).
- **(B) consensus evidence** (`cs581/bbevidence`, `cs581/protcons`): keep only the backbone residue pairs
  that both L-INS-i and FFT-NS-2 `--op 3` make.

Variants (`code/mg.py`, which builds on `gcmgen/code/gg.py` and `bbevidence/code/bbe.py`):

| name | meaning |
|---|---|
| `linsi` | control: MAGUS's own merge (reproduces MAGUS's output) |
| `linsi&fftns2-op3` | hard consensus (B) |
| `softW:linsi&fftns2-op3` | soft consensus: unconfirmed pairs keep weight W relative to confirmed ones (W = 0.5, 0.2) |
| `wsoft0.03:linsi&fftns2` | gcmgen's best general recipe when I started: soft weight 0.03, FFT-NS-2 as second opinion |
| `linsi\|cons0.7` | self-consistency mask (columns re-aligned by < 70% of the other backbones become insertion columns) |
| `sm0.5:linsi\|cons0.7` | soft mask: masked pairs keep weight 0.5 |
| `ss:X` | self-soft whose pseudo-backbones come from variant X's merged alignment (two stages) |

So `ss:linsi` is (A) as published. `ss:linsi&fftns2-op3` and `ss:softW:…` are "self-soft with
consensus-filtered evidence" (hard and soft). `ss:linsi|cons0.7` is "self-soft with a self-consistency mask".
`ss:wsoft0.03:linsi&fftns2` is the gcmgen recipe with self-soft. Each `ss:` row's runtime includes its
first stage.

## 2. Data and protocol

- **DNA/RNA (6):** ROSE 1000L1, 1000M2, 1000S1, 1000S2 (R0); RNASim 1000 R0; 16S.M R0. These are cached
  MAGUS runs from the worker branches. 1000L3 was queued but not reached.
- **Protein (9):**
  - BAliBASE RV100 BBA0039, 0067, 0101, 0134, 0154 and 0190 (cached MAGUS runs);
  - BBA0117 (fresh MAGUS draw);
  - AliSim SIMMOD_R1 and SIMHIGH_R1, regenerated exactly as in `cs581/protbench` (`code/sim.sh`; reference
    lengths match: 7348 / 11468) with fresh MAGUS draws.
  - BBA0081 was queued but not reached.
- **Protocol.** Fresh draws use `gcmx.bbtool_bench` with the paper's flags. All comparisons are paired on
  one MAGUS draw (same subsets, same backbone sequence sets). The merge flags are bbe.py's. FastSP scores
  the whole alignment.
- **Variant list was trimmed for time** (`code/variants.txt` history). Two things drove it:
  - Datasets took 20–30 min each.
  - On BBA0134 the `ss:linsi` merge took 2,292 s and `ss:linsi&fftns2-op3` was still in the trace after
    60 min.

  So I capped self-soft merges at 900 s and dropped variants that were already dead or dominated:
  - after 2 DNA sets: self-soft on hard bases, on DNA only;
  - after 3–7 sets: `soft0.2`, `ss:soft0.2` and the soft mask;
  - after 9: `ss:soft0.5`.

  Their partial rows are in the tables and are marked by n.
- **In-sample caveat.**
  - `ss:wsoft0.03:linsi&fftns2` was added after 7 datasets: 1000L1, 1000M2, RNASim, BBA0039, 0067, 0101
    and 0134 (it was then backfilled on them).
  - gcmgen chose `wsoft0.03` on overlapping BAliBASE/ROSE replicates.
  - The 8 sets run after the choice serve as the honest check: 1000S1, 1000S2, 16S.M, BBA0154, BBA0190,
    BBA0117, SIMMOD_R1 and SIMHIGH_R1.

## 3. Results

### 3.1 Per dataset (Δ error, points)

| dataset | MAGUS err | hard cons | soft0.5 | cons mask | `wsoft0.03` | **`ss:linsi`** | **`ss:wsoft0.03`** |
|---|---|---|---|---|---|---|---|
| 1000L1 | 7.47 | +7.22 | −0.01 | +9.16 | −0.04 | −0.70 | −0.62 |
| 1000M2 | 8.23 | +15.65 | +0.00 | +9.49 | +0.01 | −0.73 | −0.67 |
| 1000S1 | 9.78 | +3.89 | −0.01 | +3.77 | +0.01 | −2.30 | −2.28 |
| 1000S2 | 4.74 | +0.04 | −0.00 | +0.14 | −0.21 | −0.36 | −0.65 |
| RNASim | 9.92 | +0.58 | −0.34 | −0.08 | +0.42 | −0.31 | −0.71 |
| 16S.M | 12.86 | −0.40 | −0.07 | −0.37 | −0.64 | −0.06 | −0.34 |
| BBA0039 | 4.69 | −0.16 | −0.04 | −0.04 | −0.06 | −0.01 | +0.01 |
| BBA0067 | 26.28 | −0.96 | +0.07 | −2.47 | −0.68 | −0.31 | −1.04 |
| BBA0101 | 29.19 | −1.59 | −0.11 | −2.94 | −1.89 | −1.22 | −1.83 |
| BBA0117 | 12.61 | +0.37 | −0.00 | −0.04 | +1.07 | −0.21 | +1.02 |
| BBA0134 | 18.34 | −1.26 | +5.08 | −1.50 | −0.20 | **+11.07** (2,292 s) | −0.20 |
| BBA0154 | 21.66 | −1.76 | −0.23 | −1.05 | −1.64 | −0.45 | −1.34 |
| BBA0190 | 23.40 | −1.96 | −1.05 | −0.23 | −1.84 | +0.13 | −1.61 |
| SIMMOD_R1 | 13.62 | −1.97 | −0.19 | −1.44 | +1.33 | +1.23 | **+2.83** |
| SIMHIGH_R1 | 25.41 | −2.51 | −0.76 | +0.76 | +1.08 | timeout | timeout |

### 3.2 Per variant (SPFN / SPFP split)

| variant | group | n | Δ err | W/T/L | p | Δ SPFN | Δ SPFP |
|---|---|---|---|---|---|---|---|
| `ss:wsoft0.03:linsi&fftns2` | DNA/RNA | 6 | −0.88 | 6/0/0 | 0.031 | −0.84 | −0.92 |
| | protein | 9 | −0.27 | 5/1/3 | 0.46 | +1.94 | −2.49 |
| | pooled | 15 | −0.53 | 11/1/3 | 0.068 | +0.75 | −1.81 |
| `ss:linsi` | DNA/RNA | 6 | −0.74 | 6/0/0 | 0.031 | −0.86 | −0.63 |
| | protein | 9 | +1.28 | 4/1/4 | 0.95 | +3.12 | −0.56 |
| | pooled | 15 | +0.41 | 10/1/4 | 0.15 | +1.41 | −0.59 |
| `ss:soft0.5:linsi&fftns2-op3` (dropped at 9 sets) | DNA/RNA | 4 | −1.09 | 4/0/0 | 0.13 | −1.22 | −0.97 |
| | protein | 4 | −0.47 | 3/0/1 (BBA0134 timeout) | 0.25 | +0.08 | −1.03 |
| | pooled | 8 | −0.83 | 7/0/1 | 0.016 | −0.66 | −0.99 |
| `ss:soft0.2:linsi&fftns2-op3` (dropped) | pooled | 6 | −0.59 | 5/1/0 | 0.031 | −0.17 | −1.00 |
| `ss:linsi&fftns2-op3` (hard base) | DNA/RNA | 2 | +2.20 | 1/0/1 | – | +3.84 | +0.56 |
| | protein | 4 | −0.98 | 3/0/1 (BBA0134 timeout) | 0.25 | +0.45 | −2.41 |
| `ss:linsi\|cons0.7` (mask base) | DNA/RNA | 2 | +4.84 | 1/0/1 | – | +10.93 | −1.26 |
| | protein | 3 | −1.93 | 2/1/0 | 0.25 | +1.59 | −5.46 |
| `linsi&fftns2-op3` | DNA/RNA | 6 | +4.50 | 1/1/4 | 0.094 | +9.56 | −0.57 |
| | protein | 9 | −1.31 | 8/0/1 | 0.012 | −0.15 | −2.48 |
| `linsi\|cons0.7` | DNA/RNA | 6 | +3.69 | 2/0/4 | 0.22 | +8.58 | −1.20 |
| | protein | 9 | −0.99 | 6/2/1 | 0.027 | +0.42 | −2.41 |
| `wsoft0.03:linsi&fftns2` | DNA/RNA | 6 | −0.08 | 2/3/1 | 0.69 | +0.57 | −0.72 |
| | protein | 9 | −0.32 | 6/0/3 | 0.43 | +1.98 | −2.61 |
| `soft0.5:linsi&fftns2-op3` | pooled | 15 | +0.16 | 7/6/2 | 0.048 (worse) | +0.60 | −0.29 |
| `sm0.5:linsi\|cons0.7` (dropped) | pooled | 6 | −0.11 | 3/2/1 | 0.16 | +0.07 | −0.29 |

### 3.3 Runtime

The extra wall seconds over MAGUS are measured on this machine (4 threads): backbone realignment, masking,
both merges, the pseudo-backbones, HMM extension and the split, minus MAGUS's own merge. This is
`gcmx.e2e_bench`'s accounting: MAGUS wall plus the measured wall of the added steps.

Two caveats:
- MAGUS's own wall is measured here only for the 3 fresh draws.
- For the cached replicates it comes from the original runs on other machines, so those percentages are
  only indicative.

I did not run a separate e2e_bench pass, because the time went into the accuracy queue.

| variant | median extra s | median % of MAGUS | fresh draws (MAGUS measured here): BBA0117 / SIMMOD_R1 / SIMHIGH_R1 |
|---|---|---|---|
| one-stage consensus / mask | ≈ 0–20 | 0–3% | −2% to +13% |
| `wsoft0.03:linsi&fftns2` | 42 | 4% | +17% / +1% / +4% |
| `ss:linsi` | 169 | 23% | +26% / +25% / **> +96% (timeout at 900 s)** |
| `ss:wsoft0.03:linsi&fftns2` | 250 | 20% | +51% / +29% / **> +101% (timeout)** |

The earlier e2e pilot measured self-soft at +9% on ROSE DNA. Here the second merge is 100–230 s on
1000-sequence sets. On hard proteins it is unbounded: BBA0134 took 2,292 s, and SIMHIGH_R1 exceeded the
900 s cap for both self-soft variants. There MAGUS's minclusters trace prints "Heap limit 5000 reached ...
Setting search strategy to fully greedy".

## 4. What the data say

1. **Self-soft is the DNA/RNA ingredient, and it is robust there.**
   - Every self-soft variant on a non-catastrophic base wins on all DNA/RNA sets: 6/0/0 for `ss:linsi` and
     for `ss:wsoft0.03`, about −0.7 to −2.3 points.
   - Both SPFN and SPFP improve.
   - This replicates pilot (A).
2. **Consensus filtering is the protein ingredient, and it is unsafe on DNA**, as in (B).
   - Hard consensus and the self-consistency mask win on proteins (8/0/1, p = 0.012; 6/2/1, p = 0.027).
   - On ROSE they lose +3.7 to +15.7 points, all through SPFN: FFT-NS-2's overlap with L-INS-i is 0.02 on
     1000L1 (gcmgen's `agree` statistic).
   - The self-consistency mask fails on ROSE too (+9.2 on 1000L1/1000M2). Masking by cross-backbone
     agreement is not a DNA-safe substitute for a second aligner.
   - Self-soft on a catastrophic base only partly repairs it (1000L1: hard consensus +7.22, ss on it +4.56).
3. **Soft weighting makes consensus safe on DNA but erratic and weak on proteins.**
   - `soft0.5` is a tie on DNA (−0.07, 2/4/0) but +5.08 on BBA0134 and +0.31 overall on proteins.
   - `wsoft0.03` is −0.32 (6/0/3) on proteins.
   - The weak and erratic protein effect of soft evidence matches gcmgen's interim tables.
4. **Combining them adds up on DNA and on BAliBASE, but not on held-out proteins.**
   - `ss:wsoft0.03` keeps self-soft's DNA gain (6/0/0).
   - It gets the consensus gain on 5 of 6 non-trivial BAliBASE sets.
   - It neutralizes self-soft's BBA0134 catastrophe (−0.20 vs +11.07).
   - But on the simulated proteins self-soft itself hurts: `ss:linsi` +1.23 on SIMMOD_R1, and timeouts on
     SIMHIGH_R1. Stacking on `wsoft0.03` (itself +1.33 and +1.08 there) makes it worse: +2.83.
   - BBA0117 +1.02 comes from the base (`wsoft0.03` +1.07).
   - On proteins, self-soft raises SPFN by +3.1 points on average (`ss:linsi`). Splitting each subset into
     3 groups lets GCM under-merge, and the trace on 75 groups searches far longer. BBA0134's `ss:linsi`
     output has 7,977 columns for a 3,186-column reference.
5. **No reference-free guard I tried separates the failures.**
   - The ratio of the self-soft output length to its base output length is 1.05–1.37 on wins and losses
     alike (16S.M 1.59 is a win; SIMMOD_R1 1.16 is a loss).
   - Falling back to the base on a timeout fixes only the timeouts. `ss:wsoft0.03` would still lose on
     BBA0117, SIMMOD_R1 and SIMHIGH_R1 (via its base, +1.08).

## 5. Verdict per variant

| variant | verdict | why |
|---|---|---|
| `ss:wsoft0.03:linsi&fftns2` | **unclear, leaning not promising** | best pooled (−0.53, 11/1/3, p = 0.068); DNA 6/0/0; protein not significant, with 2 held-out losses (+1.02, +2.83) and a timeout; +20–30% runtime, unbounded on hard proteins |
| `ss:soft0.5:linsi&fftns2-op3` | unclear (n = 8) | −0.83, 7/0/1, p = 0.016 on 8 sets, but times out on BBA0134 and was not run on the held-out proteins |
| `ss:linsi` (self-soft) | not promising as a general method | DNA 6/0/0, but proteins +1.28 with a +11.07 catastrophe and a timeout |
| `ss:linsi&fftns2-op3`, `ss:linsi\|cons0.7` | not promising | inherit their base's DNA catastrophe (+4.6, +10.3 on 1000L1) |
| `linsi&fftns2-op3`, `linsi\|cons0.7` | not promising as general methods | protein winners, DNA catastrophic (known) |
| soft one-stage (`soft0.5`, `soft0.2`, `sm0.5`, `wsoft0.03`) | not promising | safe on DNA, protein gains small and erratic |

**Go/kill: kill.** No single variant wins on both data types under a Wilcoxon test without a type switch.

## 6. Suggestions to the orchestrator

- The two ingredients are complementary. Each one's failure mode is mechanistic and specific to its data
  type:
  - Consensus on DNA: the second opinion is near-random.
  - Self-soft on proteins: SPFN rises and the 75-group trace blows up.
- A general method therefore needs a reference-free regime signal rather than a better blend. Two
  candidates:
  - **The bbevidence/protcons support gate.** Here: 0.51 and 0.58 on BBA0101 and BBA0039, 0.68 and 0.78 on
    1000L1 and RNASim. It could pick the stage-1 evidence.
  - **A bounded self-soft merge.** For example, keep the original subsets as hard constraints for the
    `ss` merge on low-divergence groups, or cap the trace and fall back to the stage-1 alignment.
- Whether self-soft's protein SPFN loss is fixable (fewer groups, m = 2; or adding the original subset
  alignments as extra evidence, as `gcmx.split`'s docstring suggests) is the cheapest next experiment. It
  was out of time here.
- I did not evaluate the support gate on all 15 sets. Its τ = 0.615 was fit on overlapping BAliBASE/ROSE
  replicates.

## 7. Reproduce

```
bash cs581/code/setup.sh
for k in MOD HIGH; do bash cs581/magusgen/code/sim.sh $k 1; done
bash cs581/magusgen/code/run_all.sh      # restartable; variants re-read per job from code/variants*.txt
python3 cs581/magusgen/code/summarize.py > cs581/magusgen/results/tables.md
```

Work directories are under `/opt/work/mg` (`reps/<NAME>/variants/<variant>/out.fasta`, logs).
