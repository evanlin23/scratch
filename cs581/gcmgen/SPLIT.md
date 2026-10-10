# Train / held-out split for the adaptive switch (variant 4)

Fixed on 2026-10-10 before any agreement statistic or variant result was looked at (only MAGUS's own
control merges on BBA0101 and 1000M2 had finished).

- **Training (θ is chosen here):** BBA0101, BBA0134, BBA0067, BBA0039, SIMMOD_R1, SIMHIGH_R1 (protein);
  1000M2, 1000L1, 1000L2, 16S.M (DNA/RNA). The six phase-1 sets are all in training.
- **Held out (evaluated once with the fixed θ):** BBA0154, BBA0190, BBA0081, BBA0117, SIMMOD_R2, SIMHIGH_R2
  (protein); 1000L3, 1000M3, 1000S1, 1000S2, RNASim (DNA/RNA).

Rule: filter (use the chosen filtered/soft variant) iff overlap ≥ θ, where overlap = mean over the 10
backbones of the fraction of L-INS-i cross-subset residue pairs that the second aligner also aligns;
otherwise keep MAGUS's evidence. θ = the value maximising the mean Δ over the training sets (ties broken
by the midpoint between the adjacent training overlaps).

## Recipe selection (added 17:30, after the training results for the first variants came in)

The final "general recipe" is also selected on the training sets only: the grid
w ∈ {0.01, 0.03, 0.1} × k ∈ {1, 2, 3, 4, 5} (soft weight of unconfirmed FFT-NS-2 pairs × minimum
number of L-INS-i backbones supporting an edge) is run on the 10 training replicates; the held-out
replicates are scored only with the selected recipe and the named comparison variants. Some held-out
replicates were already run with a few grid points before this rule was written (they are reported,
but not used to choose).
