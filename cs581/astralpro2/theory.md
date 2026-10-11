# ASTRAL-Pro under pure GDL: an exact inconsistency theorem

*CS581 pilot "astralpro2", 2026-10-10. Branch `claude/cs581-astralpro2`. It builds on `../gdlcons/results/theory.md`; cand2–cand5 refer to that file.*

**Code for each result.**
- Exact values: `code/exact_proof.py` (interval arithmetic) and `code/famR.py` (quadrature).
- Simulation checks: `code/famR_sim.py` and `code/ownmap.py`.
- `code/apro.py` re-implements ASTRAL-Pro3's rooting, tagging and SQ-class counting. On 100,000 family-R trees it agrees with the ASTRAL-Pro3 binary (Section 5).

## 0. Summary

1. **Theorem 1 (exact scores).** On the 4-taxon caterpillar (((A,B)x,C)y,D), assume no events on branch x. With the correct root and ASTRAL-Pro's species-overlap tags, the limiting per-family ASTRAL-Pro score of each topology is a closed-form sum over the duplication nodes of the y-branch.
2. **Theorem 2 (sign of the bias).** At every duplication node whose two sides carry unequal numbers of y-copies, the expected contribution favours AC|BD over the true AB|CD iff b < c. Here a, b, c are the per-copy survival probabilities of A, B, C. It favours AD|BC iff a < c.
   - *Corollary 1.* If c ≤ min(a, b), overlap tagging is consistent for **every** duplication/loss process on the y-branch.
3. **Theorem 3 (rigorous counterexample).** Take a pure-birth y-branch with λT = ln 20, pure loss on the tips with a = b = 1/100 and c = 1/5, and no events on x or D.
   - The limiting scores per y-copy are correct = 2.4025·10⁻⁴ against 3.2663·10⁻⁴ for each wrong topology. These are interval enclosures to 50 digits.
   - So ASTRAL-Pro with the correct root and its own tagging rule is **statistically inconsistent** under GDL.
   - Under random hidden-paralog mislabelling it is inconsistent iff the mislabelling probability is above q* = 0.188.
4. **Proposition 4 (stock ASTRAL-Pro3, min-duplication rooting; statistical certificate, not a proof).** In the same model, the stock software's own rooting makes the deficit 2.6× larger: correct − AC|BD = −0.00487 ± 0.00004 per family (z = −117, 800,000 families). Its own rooting also creates failures where the correct root is consistent (Section 5).
5. **Proposition 5 (reconciliation is self-confirming).** Tag every gene tree by LCA reconciliation against any rooted species tree T. Then every SQ displays T, so ASTRAL-Pro returns T. Reconciliation against the true tree "fixes" the counterexample only by assuming the answer. Reconciliation against a first-pass ASTRAL-Pro tree keeps the first-pass tree, including its error.

## 1. Model and definitions

**Species tree.** S = (((A,B)x,C)y,D), rooted.
- A gene family starts as one copy at the root of S. There is no root branch.
- On every branch e, each copy duplicates at rate λ_e and is lost at rate μ_e (linear birth–death); this is pure GDL.
- At a speciation each copy enters both child branches.
- The gene tree G is the locus tree with lost lineages pruned and degree-2 nodes suppressed. Its leaves are labelled by species.

**Branch notation.**
- y-branch: the branch from the root to y, with length T and rates (λ, μ).
- *y-copies*: the copies that reach node y.
- For X ∈ {A, B, C}, **s_X** is the probability that a single copy entering the branch above X leaves at least one descendant at X. We write a = s_A, b = s_B, c = s_C and x̄ = 1 − x.

**ASTRAL-Pro** (Zhang et al. 2020, MBE 37:3292; definitions as quoted in `../gdlcons/results/prior_art.md`):
- **Tagging.** A node can be tagged S only if its children's species sets are disjoint. The *species-overlap* (ovl) tagging tags S exactly those nodes.
- **SQ.** Four leaves from four species such that the LCA of every three of them is tagged S.
- **Anchor LCA.** The LCA of the two degree-3 nodes of G|Q. Two SQs are equivalent iff they have the same anchor LCA.
- **Score.** The per-family score of a quartet topology is the number of equivalence classes displaying it. On 4 taxa, ASTRAL-Pro returns the topology with the largest total score.
- **ASTRAL-Pro3's rooting** (ASTER `src/astral-pro.cpp`, `scoreSubtree` and `annotateTree`):
  - every edge is tried as the root;
  - the root minimising Σ_{D-nodes v} (1 + [L_v ≠ L_v ∪ R_v] + [R_v ≠ L_v ∪ R_v]) is kept, where L_v and R_v are the children's species sets;
  - ties are broken uniformly at random.

**Consistency.** ASTRAL-Pro is consistent if P(output = S) → 1 as the number of i.i.d. families N → ∞.
- Per-family scores are bounded by the number of internal nodes, which has finite mean in a linear birth–death process on a finite tree.
- So by the strong law of large numbers, for 4 taxa: if some wrong topology has a strictly larger *expected* per-family score than AB|CD, then P(output = AB|CD) → 0.
- If AB|CD has the strictly largest expected score, then P(output = AB|CD) → 1.
- Families with fewer than four species score 0 for every topology, so filtering them out does not change the comparison.

## 2. Where classes live (correct root)

**Lemma 1.** Root G correctly and tag it with any tagging in which the root is S. Then the SQ classes are in one-to-one correspondence with the nodes w in the ABC part of G such that:
- w is tagged S;
- the species sets of w's two children, restricted to {A, B, C}, form a 2|1 partition {P, {A,B,C}∖P} of {A,B,C}.

The class at w displays P | ({A,B,C}∖P) ∪ {D}.

*Proof.*
- The root of G is the root speciation. Its children are the D part and the ABC part, so every LCA of a triple that contains a D-leaf is the root, which is S.
- Hence Q = (a, b, c, d) is an SQ iff w := LCA(a, b, c) is tagged S.
- G|Q is the caterpillar (((·,·),·),d) whose cherry lies inside one child of w. Its degree-3 nodes are the cherry's LCA and w, so the anchor LCA is w.
- A class therefore exists at w iff some choice of a, b, c below w has its three LCAs at w with all three triples' LCAs S. That is exactly the 2|1 condition on w's children's species sets, plus w tagged S.
- Under ovl tags, an S node has disjoint children. So all SQs anchored at w put the same pair P on one side, and the class has one topology. ∎

## 3. Model R: exact scores and the sign of the bias

**Model R.**
- Branch x has no events.
- Branches A, B, C and D are arbitrary; they enter only through s_A, s_B, s_C and P(D survives) > 0.
- The y-branch is an arbitrary birth–death process.

**Lemma 2 (factorisation).** In model R, condition on the y-branch duplication tree, i.e. on how the y-copies are joined by duplications. Then:
- the species sets of the y-copies (restricted to {A, B, C}) are i.i.d.;
- inside each copy, A, B and C are present independently with probabilities a, b and c.

So a subtree that contains j ≥ 1 y-copies has species set S with probability

  p_j(S) = Π_{X∈S} (1 − x̄^j) · Π_{X∈{A,B,C}∖S} x̄^j.

*Proof.* With no events on x, a y-copy sends exactly one copy into each of the A, B and C branches. Those three branch processes are independent birth–death processes, and they are independent across y-copies. ∎

**Which nodes qualify under ovl tags.**
- *y-copies.* A y-copy has children with species sets ⊆ {A,B} and ⊆ {C}, so it is always S. It carries the AB|CD class iff all of A, B and C survive. Expected number: E[N_y]·abc.
- *Nodes in x and the tips.* Their sets lie inside {A,B} or inside one species, so they carry no class.
- *Duplication nodes w on the y-branch.* Let i and j be the numbers of y-copies on w's two sides (if either is 0, w is suppressed). Under ovl tags, w carries a class with pair P iff its two sides have exact species sets P and {A,B,C}∖P. These are the "hidden paralogs".

**Theorem 1 (limiting scores, ovl tags, correct root).** In model R, writing π_D = P(D present) and V_P for the pair-P term (V_AB, V_AC, V_BC):

- E[score_{AB|CD}] = π_D · ( E[N_y]·abc + E Σ_w V_AB(i_w, j_w) )
- E[score_{AC|BD}] = π_D · E Σ_w V_AC(i_w, j_w)
- E[score_{AD|BC}] = π_D · E Σ_w V_BC(i_w, j_w)

where V_P(i, j) = p_i(P)·p_j(P^c) + p_j(P)·p_i(P^c) and P^c = {A,B,C}∖P. With true tags the two wrong scores are 0, which is Zhang et al.'s Theorem 2.

*Proof.* Combine Lemma 1, Lemma 2 and the node list above. D's presence is independent of the ABC part, which gives the factor π_D. The sum over nodes is a sum of indicators, so linearity of expectation applies, conditioning on the duplication tree. ∎

**Theorem 2 (sign of the per-node bias).** Let 1 ≤ i < j and 0 < a, b, c < 1. Then

  V_AB(i, j) − V_AC(i, j) = (ā^j − ā^i) · g_B · g_C · (r_B − r_C),

where, for X ∈ {B, C}, with u = x̄^i and v = x̄^j:
- g_X = u(1 − v) > 0;
- r_X = (1 − u)v / (u(1 − v)) = (y^m − y^j)/(1 − y^j), with y = x̄ and m = j − i.

The map y ↦ r is strictly increasing on (0, 1). Hence **sign(V_AB − V_AC) = sign(b − c)**. In the same way, sign(V_AB − V_BC) = sign(a − c). For i = j both differences are 0.

*Proof.*
1. Set f_X = (1 − u_X)v_X and g_X = u_X(1 − v_X). Then p_i(AB)p_j(C) = f_A f_B g_C, and the other products follow the same pattern. This gives V_AB = f_A f_B g_C + g_A g_B f_C and V_AC = f_A g_B f_C + g_A f_B g_C.
2. Subtracting, V_AB − V_AC = (f_A − g_A)(f_B g_C − f_C g_B), with f_A − g_A = v_A − u_A = ā^j − ā^i < 0.
3. Factor f_B g_C − f_C g_B = g_B g_C (r_B − r_C).
4. Monotonicity: d r/dy has the sign of y^{m−1} φ(y), where φ(y) = m − j y^i + i y^j. Now φ(1) = 0 and φ′(y) = ij y^{i−1}(y^m − 1) < 0, so φ > 0 on (0, 1).
5. So r_B > r_C ⇔ b̄ > c̄ ⇔ b < c, and in that case V_AB − V_AC < 0. ∎

**Interpretation.** At a duplication whose two sides carry unequal numbers of copies, the side with more copies is more likely to contain *every* species, and that inflates whichever pattern pairs the best-surviving outer taxon C with a cherry member. The bias pulls the worse-surviving of A and B toward the outgroup D.

**Corollary 1 (sufficient condition).** In model R, if c ≤ min(a, b) and abc > 0, then ASTRAL-Pro with the correct root and ovl tags is consistent, whatever the rates and length of the y-branch. Every hidden-paralog term then favours AB|CD and O = π_D E[N_y] abc > 0. The same holds under random hidden-paralog mislabelling `d2sv(q)`: scale every V by q.

**Closed form for a pure-birth y-branch** (rate λ, length T, μ = 0).
- A duplication at time t has two independent subtrees, each with Geometric(q) many y-copies, where q = e^{−λ(T−t)}.
- The density of duplications at time t is λe^{λt}. Substituting q gives

  E[score_P]/(π_D e^{λT}) = [P = AB]·abc + 2 ∫_{q0}^{1} p_q(P) p_q(P^c) dq,  q0 = e^{−λT},

  where p_q(S) = Σ_{U⊆S} (−1)^{|S∖U|} G(z_U), with z_U = Π_{X∉U} x̄ and G(z) = E[z^k] = q/(q + (1 − z)/z).
- The integrand is rational in q, so the integral is an explicit combination of rationals and logarithms (`exact_proof.py` gives it term by term).

**Theorem 3 (counterexample).** Take λT = ln 20 (q0 = 1/20), a = b = 1/100, c = 1/5, and no events on x and D. For example: tips with pure loss, μ_A t_A = μ_B t_B = ln 100 and μ_C t_C = ln 5. Interval arithmetic at 50 digits (`python exact_proof.py 1/100 1/100 1/5 1/20`) gives, per y-copy:

| quantity | value |
|---|---|
| O = abc | 2.0000·10⁻⁵ |
| H_AB | 2.20250·10⁻⁴ |
| H_AC = H_BC | 3.26632·10⁻⁴ |
| score(AB\|CD) − score(AC\|BD) | **−8.6382·10⁻⁵** (interval width < 10⁻⁵⁰) |

So the true topology has the smallest of the three limiting scores, and **ASTRAL-Pro with the correct root and its own overlap tags is statistically inconsistent**.
- For random hidden-paralog mislabelling at rate q, it is inconsistent iff q > q* = O/(H_AC − H_AB) = 0.188.
- *Simulation check* (`famR_sim.py`, 800,000 families at λT = 3): per family, AB|CD = 0.00478, AC|BD = 0.00663 and AD|BC = 0.00647. The formula predicts 0.00483, 0.00657 and 0.00657.

**Phase region.** Results from `famR_minlam.py` and `famR_phase.py` (exact formula, a = b ≥ 10⁻³):
- Failure needs λT ≳ 1.4, i.e. about 4 or more expected y-copies. The minimum relative margin is +0.034 at λT = 1, −0.017 at λT = 1.5, −0.059 at λT = 2 and −0.118 at λT = 3.
- It needs a, b ≲ 0.05: with a = b = 0.1 there is no failing c at any λT.
- C must survive much better, but not too well: for λT = 3, c ∈ (0.027, 0.51) when a = b = 0.01, and c ∈ (0.08, 0.34) when a = b = 0.05.

**Relation to the earlier pilot.** cand2 and cand5 had events on x and on tips, so they lie outside model R. Their failures came from the same mechanism and matched the general formula numerically. Model R isolates it and gives a closed form plus an exact sign law.

## 4. Proposition 5: reconciliation-based tagging is self-confirming

**Proposition 5.** Let T be any rooted binary species tree on the gene tree's species. Root G anywhere, and tag each node v as D iff its LCA image M(v) equals M(c) for a child c. Then every SQ of G displays T restricted to its four species. Consequently, for 4 taxa, ASTRAL-Pro returns T whenever at least one SQ exists.

*Proof.*
1. Let u = LCA(x, y, z) be tagged S, with x and y below one child u₁ of u and z below the other child u₂.
2. Since u is S, M(u₁) ≠ M(u) ≠ M(u₂). If M(u₁) and M(u₂) were in the same child subtree of M(u), their LCA would lie strictly below M(u), a contradiction. So they lie in different children.
3. Hence LCA_T(s_x, s_y) ≤ M(u₁) < M(u) = LCA_T(s_x, s_z), so T displays the rooted triplet s_x s_y | s_z.
4. In an SQ all four triples have S LCAs, so the rooted restriction G|Q agrees with T on all four triplets. A rooted 4-leaf tree is determined by its triplets, so G|Q = T|Q. ∎

**Consequences.**
- "Reroot and retag by reconciliation against a first-pass species tree" makes ASTRAL-Pro output the first-pass tree. Every tree is a fixed point of the iteration, so it cannot correct a wrong first pass.
- Against the true tree it is trivially consistent, but that is an oracle, not a method. In model R, true-tree reconciliation tags the AC/B and BC/A hidden paralogs D and keeps the AB/C ones S (checked in `runmeth.py`, method `recon_true`).
- The bias in Theorems 1–3 sits in which quartets survive a topology-blind filter. Any fix must either use information beyond the leaf-labelled topology (branch lengths and dating, sequence data) or avoid filtering, as ASTRAL-multi does. Whether ASTRAL-multi survives model R is tested empirically; Legried et al.'s proof assumes [M] constant rates.

## 5. Stock ASTRAL-Pro3 rooting (min-duplication/loss score)

I found no tractable closed form for the expectation under min-score rooting: the chosen root depends on the whole tree, and ties are broken at random. What follows is a numerical certificate, not a proof.

**Validation of `apro.py`.**
- Over 100,000 family-R trees, totals of tie-averaged own-rooting scores (AB/AC/AD): ASTRAL-Pro3 binary 944 / 1,546 / 1,518; `apro.py` 965 / 1,523 / 1,520. The differences are within the tie-break noise.
- In the Theorem 3 configuration the true root is among the optimal roots in 97.0% of families, and 30.6% of families have tied optima.

**Theorem 3 configuration, 4 × 200,000 families.** Per family (between-seed SE):

| rooting/tags | AB\|CD | AC\|BD | AD\|BC | AB − AC |
|---|---|---|---|---|
| correct root + ovl | 0.00478 | 0.00663 | 0.00647 | −0.00184 ± 0.00005 |
| **stock own rooting** | 0.00817 | 0.01304 | 0.01297 | **−0.00487 ± 0.00004** |

**Own rooting enlarges the failure region** (`ownmap.py`, family R, 100,000 families per cell; relative margin = (correct − best wrong)/total):

| a = b | c | λT | exact, correct root | stock own rooting (z) |
|---|---|---|---|---|
| 0.02 | 0.4 | 2 | +0.003 | **−0.089 (z = −9.0)** |
| 0.02 | 0.6 | 2 | +0.129 | **−0.151 (z = −8.7)** |
| 0.05 | 0.6 | 1 | +0.295 | +0.076 |
| 0.1 | 0.6 | 1 | +0.356 | +0.137 |

The full table is in `results/ownmap.jsonl`.

**The two wrong-tree mechanisms are different.**
- Correct-root overlap tagging fails only when C survives better, but not too much better (Theorem 2 together with a bounded c-range).
- Own rooting fails increasingly as c grows. Min-score rooting puts the root inside the abundant C lineage, as cand4 showed in the earlier pilot. Then hidden-paralog-like S nodes appear on the re-rooted path, and they favour pairing C with the cherry.

**Open (and the natural theorem for a project).**
- An exact formula, or a provable bound, for the min-score rooting in model R.
- Two ingredients are already in hand: (i) the true root is optimal in about 97% of families; (ii) by Claim 1 of Zhang et al., rerooting along an all-S path does not change the score.
