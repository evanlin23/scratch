# ASTRID-Pro, round 2: theorem, broader benchmark, scaling

*CS581 pilot, 2026-10-10. Branch `claude/cs581-astridpro2`. Code is in `code/`, raw results in `results/`, the
theorem in `theory.md`, and the pre-registration in `PREREG.md`. Machine: 4 cores and 15 GB RAM. All methods ran
single-threaded, but up to 4 jobs ran at once (see the runtime caveat).*

__VERDICT__

## 1. Theory (full statement and proofs in `theory.md`)

**Setting.** Pure GDL (any branch-specific, time-varying birth–death rates; no ILS), true gene trees, true root and
D/S tags. The ASTRID-Pro distance counts the speciation nodes on the path between orthologous copies, including the
LCA and the gene-tree root.

- **Thm 3 (exact limit).** The averaged matrix converges a.s. to a metric that is *additive on the species tree S*.
  Internal edge c has length

  β(c) = ½ [ s(c₁) + s(c₂) + s(sib c) − s(c) ]

  where s(·) is the survival probability of a copy entering a branch. These probabilities have closed forms for
  linear birth–death. Leaf edges are ½[1 + s(sib) − s(leaf)] ≥ 0.
- **Thm 4 (iff).**
  - If β(c) > 0 at every interior node, ASTRID-Pro + NJ/FastME is statistically consistent.
  - If some β(c) < 0, the limit violates the four-point condition. The true split around c becomes the *largest*
    quartet sum, so every distance method fails on that quartet.
  - This corrects the previous note, which said the limit is "always a tree metric".
- **Thm 5 (sufficient).** β(c) > 0 whenever the branch above c is not supercritical (λ ≤ μ). Leaf branches and root
  branches can have any rates. With no loss, β ≡ 1.
- **Counterexample.** Exact β = −0.267. Simulation agrees with the predicted matrix to 0.003 at 100,000 families,
  and ASTRID-Pro converges to the wrong tree. A second, random 6-taxon instance with β = −0.034 is also wrong at
  20,000 families.
- **Benchmarks are covered.** The exact minimum interior β on every DISCO and FastMulRFS simulation condition is
  ≥ 0.46, including the supercritical loss/dup = 0.5 conditions. Over 20,000 random-rate configurations on 4–8
  taxa, 1.15% fail, and all of those have a supercritical branch.
- **Prop 6 (correction).** Reweighting each counted node by 1/ŝ(off-path clade), with ŝ estimated by reconciling to
  a correct rooted first-pass tree, gives a limit with unit edge lengths. That is consistent under *any* GDL rates.
  It fixes both failing instances given a correct first pass.
- **Open.** Estimated tags (hidden paralogy), MinDup rooting, ILS, and ASTRID-multi itself. Counting the gene-tree
  root is required by the theorem with true roots. It is harmless or slightly harmful on estimated trees; see
  `astrid-pro-r0`.

__RESULTS__
