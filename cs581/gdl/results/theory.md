# Theory note: is the speciation-only ortholog distance additive under pure GDL?

## Setting
- **Species tree.** S is rooted and binary.
- **Gene families.** Each family evolves by a birth–death process on copies. Branch e has its own duplication rate λ_e and loss rate μ_e. At a species node, every copy present splits into one copy in each child branch.
- **No ILS.** The gene-family tree equals the locus tree. Extinct lineages are pruned and unifurcations suppressed.
- **Ideal input.** Gene trees are true and correctly rooted, and their nodes carry the true duplication/speciation (D/S) tags.

**ASTRID-Pro distance.** Take a pair of copies (a, b) from species A ≠ B. The pair is *orthologous* if its LCA in the gene tree is an S node, which is then the image of L = LCA_S(A, B). Define d(a, b) as the number of S nodes on the gene-tree path, counting the LCA. The per-gene distance δ̂(A, B) is the mean of d over orthologous pairs. The species distance is the mean of δ̂ over the genes that contain such pairs. **The root must be counted.** It is a speciation node at the species root. If it is not counted (ASTRID's unrooted convention), every pair crossing the root split loses 1. That shortens the root-split edge of the limiting tree metric and can make it negative. Example: a balanced quartet ((A:1,B:0.1),(C:1,D:1.5)) with GDL rates. Without the root: AB|CD = 1.87 > AC|BD = 1.76, wrong. With the root: AB|CD = 2.00 < 3.74/3.76, correct. The code exposes this as `count_root=True` (variant `pro-mean-truetags-countroot` in `probe.py`).

## Lemma (expected distance)
For an internal species node v strictly inside path_S(A, B) other than L, let off(v) be the child of v that is *not* on the path. Define s(c) as the probability that a single copy entering the branch above c leaves at least one descendant at the leaves.

The copy lineage leading to a passes through v and spawns a daughter copy into off(v). The node shows up as an S node on the observed path iff that daughter survives. That event depends only on the process inside the subtree below off(v). This subtree is disjoint from the subtrees that decide whether a and b exist. Hence, for every orthologous pair,

  E[d(a, b) | (a, b) orthologous] = 1 + Σ_{v ∈ path_S(A,B), v internal, v ≠ L} s(off(v)) =: δ(A, B).

Consequences:
- Duplications on the path never contribute, because they are not counted.
- The per-gene mean has the same expectation, so the averaged matrix converges to δ.
- This holds only for families that are **not conditioned** on how many species they contain. In practice families with fewer than 4 species are dropped, which breaks the independence slightly; see the caveats.

## Four-point analysis
Take any four leaves. They induce either a balanced or a caterpillar rooted quartet.

For **segment sums**, P_* denotes the sum of s(off(v)) over the internal nodes of one segment. s_xA is the survival of a copy entering x's child branch toward A, and s_yC is defined the same way toward C. s_yx is the survival of a copy entering the branch from y toward x.

- **Balanced, ((A,B)x,(C,D)y)r.**
  - AB|CD = P_A + P_B + P_C + P_D + 2.
  - AC|BD = AD|BC = (same) + s_xA + s_xB + s_yC + s_yD + 2P_xr + 2P_yr.
  - So the correct split is strictly smallest and the other two sums are equal.
- **Caterpillar, (((A,B)x,C)y,D).**
  - AB|CD = P + 2 + s_yx.
  - AC|BD = AD|BC = P + 2 + 2P_xy + s_xA + s_xB + s_yC.
  - Here P collects the common terms.

**In every case the two largest sums of the alternatives coincide.** So δ always satisfies the four-point condition (Buneman), and δ is **always a tree metric**. Its topology is the species tree iff every induced caterpillar satisfies

  (★)  s_yx < 2P_xy + s_xA + s_xB + s_yC.

**Sufficient condition: no supercritical branch, i.e. λ_e ≤ μ_e for every e.**
- s_yx needs at least one copy to reach some subtree that hangs off the y→x path, or one of x's two child subtrees.
- A union bound gives s_yx ≤ Σ_w m_w·s(w), where m_w is the expected number of copies reaching w.
- With λ_e ≤ μ_e on every branch, m_w ≤ 1. Then s_yx ≤ P_xy + s_xA + s_xB, which is strictly smaller than the right side of (★) because s_yC > 0.
- Hence, if s > 0 everywhere, δ is additive on the true species tree with positive internal edge lengths.
- NJ and FastME are consistent on additive matrices, so **ASTRID-Pro with ideal rooting and tagging is statistically consistent under pure GDL whenever no branch is supercritical**. This includes the common simulation setting λ = μ, as in the FastMulRFS data and the DISCO default.

**Counterexample: a supercritical branch.**
- Take (((A:1,B:1)x:1,C:1)y:1,D:3).
- On the branch above x, set λ = 3, μ = 0. On A, B and C set λ = 0, μ = 3. All other branches have no events.
- Then s_xA = s_xB = s_yC = e⁻³ ≈ 0.05. About e³ ≈ 20 copies reach x, so s_yx ≈ 0.7. Hence (★) fails: 0.7 > 0.15.
- Prediction: AC|BD = AD|BC = 1 + 3e⁻³ = 1.149 and AB|CD = 1 + s_yx. Here d_CD = s_yx because the gene-tree root is not counted.
- **Simulation** (`results/counterexample.jsonl`; true gene trees, true tags; families with at least 2 species). The numbers below are without root counting; with root counting every sum is exactly 1 larger (2.695 vs 2.1515/2.1508, `results/counterexample_countroot.jsonl`), so the conclusion is the same:
  - 20,000 families: AC|BD = 1.1515, AD|BC = 1.1508, AB|CD = 1.695.
  - The wrong split is the strict minimum and the other two sums are equal, as predicted. FastME returns the wrong tree (FN = 1) at both 2,000 and 20,000 families.
- On the same families these methods are all correct (FN = 0): ASTRAL-Pro (Zhang et al. 2020 prove it must be), ASTRID-multi, ASTRID-DISCO, and ASTRID-Pro with *inferred* tags.
- **So the ideal speciation-only ortholog distance is NOT consistent under unrestricted GDL with rate variation across branches.** It is consistent when no branch is supercritical.

## Caveats and open items
1. **Conditioning.** Real pipelines keep only families with at least 4 species. That conditions on survival in the off-path subtrees, so the Lemma becomes approximate.
2. **Estimated tags.** Species-overlap tagging calls a "hidden" duplication an S node when its two children have disjoint species sets. In the counterexample above, the version with inferred tags happened to be correct. That is luck, not a guarantee.
3. **Closest-copy (min) variant and ASTRID-multi.** The argument does not cover them. ASTRID-multi averages over paralogous pairs whose LCA is a duplication above L, and the expectations then depend on copy numbers. See the search results in `results/search_summary.md`.
4. **ILS (DLCoal).** Not covered.
5. **A possible correction.** Under (★)-violating rates the metric is still a tree metric, but with a wrong topology. That suggests a reweighting by estimated s(·) (from copy-number and presence data) could restore additivity: δ(A,B) − 1 is linear in the unknown s values. This is untested.
