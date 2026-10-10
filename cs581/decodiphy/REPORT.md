# Identifiability and choosing k in phylogenetic distance deconvolution (DecoDiPhy) — overnight pilot

*CS581 project exploration, Oct 2026. Code: `code/`, result tables: `results/`.*

## 1. Question and source open problem

Arasti, Şapcı, Rachtman, El-Kebir & Mirarab, *Deconvolving Phylogenetic Distance Mixtures*, RECOMB 2026
(bioRxiv 10.64898/2026.01.18.700179, PMC12871782; code github.com/shayesteh99/DecoDiPhy, data
github.com/shayesteh99/DecoDiPhy-Data) define **phylogenetic distance deconvolution (PDD)**: given a reference
tree R on n leaves and a mixture distance vector d (d_r = Σ_q p_q d_F(r, q)), find k query placements
(edge e_q, relative position x_q ∈ (0,1), abundance p_q, mean pendant length ȳ) with d = D·A(e)·p + C·L·A(e)·w + ȳ·1.

They leave open, verbatim:

> "Theoretical questions also remain, including characterizing unidentifiable cases, studying tolerance to
> error in d̂ …, and understanding optimality. We used a simple criterion to choose k, leaving more
> sophisticated stopping criteria for future work."

and on identifiability:

> "Claim 3. Even for fixed (true) placement edges, the relative positions x_i and proportions p_i are
> unidentifiable under some corner conditions for queries placed on adjacent edges." …
> "Even with Assumption 1, the problem is still not always identifiable … we suspect (with no proof) that
> conditions that break identifiability require extreme cases of symmetry similar to Figure 1d that are not
> expected on real data. … We leave a full characterization of identifiability to future work."

Their k rule: grow k greedily and stop at the first k with some p_q < min_p (default 0.01 for distance
input, 1000/|X| for read input), returning k−1.

Pilot questions: (1) characterize unidentifiable inputs exactly (small trees, exhaustive); (2) prove
something for k = 2; (3) test BIC / elbow / theory-motivated stopping rules against the paper's rule on
the paper's simulation protocol.

## 2. Prior art / novelty check (searched 2026-10-10)

Searched: "phylogenetic distance deconvolution" identifiability; DecoDiPhy follow-ups / citations;
DecoDiPhy stopping criterion; MISA (Balaban & Mirarab 2020, *Phylogenetic double placement of mixed
samples*, Bioinformatics 36:i335 — the k=2 special case with a k-mer distance; no identifiability
characterization found); tree distance/Laplacian identities; phylogenetic Kantorovich–Rubinstein / EMD
on trees (Evans & Matsen 2012); identifiability from average distances on networks (Xu & Ané 2023,
J Math Biol, arXiv 2110.11814 — analogous flavour: parameters near 4-cycles are not identifiable).

Found: no paper citing or extending the PDD work (it appeared Jan 2026, RECOMB May 2026); no
characterization of PDD identifiability and no alternative stopping rule. The algebraic tool we use is
the classical tree identity L·D + 2I = (2·1 − deg)·1ᵀ (Graham–Lovász for unweighted trees; weighted
version e.g. Bapat, Kirkland & Neumann 2005, LAA 401; Goubko & Veremyev 2021 use it as a constraint).
Its application to PDD appears new. **Conclusion: the open problem is still open as of today.**
(The authors' repository contains `check_identifiability`, which flags placements within 2 edges of each
other, and their simulator avoids drawing a query whose parent or grandparent is already taken. So they
are aware adjacency matters, but the paper treats it as a corner case.)

## 3. Main theoretical results

Notation: R binary, unrooted, positive branch lengths; V its nodes, I its n−2 internal nodes;
F ∈ R^{n×V}, F_{r,v} = d_R(r,v); L the weighted Laplacian of R with conductances 1/l_e.
For a placement set S (k distinct edges) write V(S) for the set of endpoints.

**Lemma 1 (barycentric reduction).** A query at relative position x on edge e=(a,b) (x measured from a)
contributes p·[(1−x)·F_{·,a} + x·F_{·,b}] to d (distance to any leaf is linear along an edge). Hence
d = F·m + ȳ·1, where the **node measure** m = Σ_q p_q[(1−x_q)δ_{a_q} + x_q δ_{b_q}] is ≥ 0, sums to 1, and
supp(m) = V(S). *All the information in d is (m, ȳ).*

**Lemma 2 (kernel).** The map Φ(m, ȳ) = (F m + ȳ 1, 1ᵀm) from R^V × R to R^{n+1} is onto, and
ker Φ = span{(−Lδ_v, −1) : v ∈ I} (dimension n−2).
*Proof.* (i) For internal v with neighbours u_1,u_2,u_3, a leaf r in branch j satisfies
Σ_i (d(r,u_i) − d(r,v))/l_i = −1 + 1 + 1 = 1, so F(−Lδ_v) = 1, and 1ᵀLδ_v = 0. (Equivalently
L D + 2I = (2·1−deg)1ᵀ.) These n−2 vectors are independent because the grounded Laplacian L[I,I] is
nonsingular. (ii) Onto: if (c, s) annihilates the image then φ(z) = Σ_r c_r d(r,z) is constant (= −s) on R;
along the pendant edge of leaf r its slope is 2c_r − Σc = 2c_r, so c = 0 and s = 0. Dimension count:
(|V|+1) − (n+1) = n−2. ∎ (Also verified numerically, `code/check_kernel.py`, 180 random trees n=4..12.)

**Theorem 1 (characterization).** Two PDD solutions produce the same d **iff** their node measures differ
by a Laplacian move: m' = m − Lα, ȳ' = ȳ − Σ_v α_v for some α supported on internal nodes. A node
measure m is realizable by edge set S' iff supp(m) = V(S') (split each node's mass among its S'-edges);
the (p,x) fibre over a fixed m and S' has dimension 2|S'| − |V(S')| (= number of shared endpoints).

Consequences (all proved by Theorem 1):

* **C1 (adjacency is never identifiable; sharpens Claim 3).** If two placements share a node, then for
  *every* parameter value the (p, x) on the true edges form a continuum of dimension ≥ 1: the pair's
  contribution depends on p_1, p_2 only through p_1+p_2 (for edges (v,a),(v,b): leaves beyond a see
  −w_1+w_2, beyond b see w_1−w_2, the rest w_1+w_2, with w = p·t measured from v). Not a corner case.
* **C2 (k is never identifiable from above).** If ȳ > 0, then for small ε > 0 the move α = ε δ_v at any
  internal v ∈ V(S) adds mass on all three neighbours of v: d is reproduced *exactly* by S plus the other
  edges of the "star" at v (k+1 or k+2 placements), with ȳ' = ȳ − ε. So zero residual persists above the
  true k; in noisy data the cheapest spurious placements are star splits next to existing placements.
  Only the **minimal** k with an exact fit is meaningful.
* **C3 (k = 1 is identifiable)** against any alternative with k' ≤ 1. Proof: with g = m' − m = −Lα, the
  negative entries of g lie in V(S) and the positive ones in V(S'), each an adjacent pair. Every leaf of the
  subtree spanned by supp α has two private outside neighbours where g has the same nonzero sign
  (objects defined in Theorem 2). They are siblings, so not adjacent, and cannot both lie in one adjacent pair.

**Theorem 2 (k = 2).** If the two placements are not adjacent (S a matching), then (S, p, x, ȳ) is
identifiable: (a) no continuum, and (b) no alternative with k' ≤ 2.
*Proof of (a).* Suppose α ≠ 0 with supp(Lα) ⊆ V(S) = {a_1,b_1,a_2,b_2}. Let U = supp α and T_U the smallest
subtree containing U. A leaf v of T_U has ≥ 2 neighbours outside T_U that are adjacent to no other U-node,
so g = −Lα equals α(v)/l ≠ 0 there; they are siblings (common neighbour v), hence non-adjacent, hence in
different S-edges (say a_1, a_2). If |U| = 1, g ≠ 0 on v and all three neighbours ⇒ 4 nodes of V(S)
including v, its S-partner and two siblings that must be S-partners of each other — impossible (siblings
are not adjacent). If |U| ≥ 2, a second T_U-leaf v' gives siblings in V(S) disjoint from the first pair,
i.e. b_1, b_2, and then v–a_1–b_1–v'–b_2–a_2–v is a cycle in a tree. ∎
*Proof sketch of (b).* Same objects; negative entries of g lie in V(S), positive in V(S'); each T_U-leaf
gives a same-sign sibling pair. |U| = 1: three siblings of one sign would need 3 edges. Two pairs of the
same sign give the 6-cycle above. The remaining case (two T_U-leaves with opposite signs, T_U a path) we
have not closed by hand; it is covered by the exhaustive check below (no counterexample, all unlabeled shapes n ≤ 12).

**Theorem 3 (closed claw, k ≥ 3): claw ⇒ non-identifiable is proved; the converse at k = 3 is verified
exhaustively (n ≤ 9).** Call v ∈ I a *closed claw* of S
if v ∈ V(S) and all three neighbours of v are in V(S) (needs ≥ 3 placements, all pairwise non-adjacent
possible). If S has a closed claw then (i) the true (p,x,ȳ) lie in a continuum and (ii) on an open set
of parameters there is a different edge set S' with |S'| = k fitting d exactly.
*Proof.* α = εδ_v moves mass only inside {v} ∪ N(v) ⊆ V(S): continuum for small |ε|. For (ii) take ε < 0
(ȳ grows, always allowed): v gains, its neighbours lose mass at rates 1/l_j; if v's S-partner u_1 has the
smallest m(u_j)·l_j, it hits 0 first and S' = S − (v,u_1) + (v,u_2) realizes the new measure. ∎
This kills the paper's "extreme symmetry" conjecture: the counterexamples form an open set of
parameters, with generic branch lengths, as soon as k ≥ 3.

**General criterion (conjecture, exhaustively verified).** For S a matching the following are
equivalent: (A) some parameters admit an exact alternative with k' ≤ k; (B) the true solution has a
continuum; (C) ker L[V ∖ V(S), I] ≠ {0}, a one-line rank test. (B) ⇔ (C) is immediate from Theorem 1 (a continuum on
fixed edges is exactly a nonzero α whose Laplacian vanishes off V(S)), so the open part is (A) ⇔ (B).
For k ≤ 3, (C) ⇔ closed claw. At k = 4 a second pattern appears, the *split claw* (§4a), but (A) ⇔ (B)
still holds in every case enumerated.

## 4. Exhaustive small-case evidence

### 4a. Configuration level (parameter-free; all unlabeled shapes; all matchings S, all S' with |S'| ≤ |S|)

`code/config_enum.py`: by Theorem 1, "S' is an exact alternative for some (p,x,ȳ)" is the feasibility of
a sign pattern on g = −Lα (zero off V(S)∪V(S'), > 0 on V(S')∖V(S), < 0 on V(S)∖V(S')), an LP; checked
for unit and random branch lengths. Full table: `results/config_enum.md`.

Aggregated over all unlabeled shapes (`results/config_enum.md`, `results/config_enum_k2_n10-12.md`):

| n | k | shapes | matchings S | S with continuum | S with exact alternative (unit / random lengths) | alternative but no closed claw | closed claw | (A)⇔(B) violations |
|---|---|---|---|---|---|---|---|---|
| 4–12 | 1 | all 82 | 1,500 | 0 | 0 / 0 | 0 | 0 | 0 |
| 4–12 | 2 | all 82 | 11,372 | 0 | 0 / 0 | 0 | 0 | 0 |
| 5 | 3 | 1 | 4 | 4 | 4 / 4 | 0 | 4 | 0 |
| 6 | 3 | 2 | 40 | 20 | 20 / 20 | 0 | 20 | 0 |
| 7 | 3 | 2 | 112 | 28 | 28 / 28 | 0 | 28 | 0 |
| 8 | 3 | 4 | 480 | 80 | 80 / 80 | 0 | 80 | 0 |
| 9 | 3 | 6 | 1,320 | 148 | 148 / 148 | 0 | 148 | 0 |
| 6 | 4 | 2 | 4 | 4 | 4 / 4 | 0 | 4 | 0 |
| 7 | 4 | 2 | 52 | 52 | 52 / 52 | 4 | 40 | 0 |
| 8 | 4 | 4 | 400 | 292 | 292 / 292 | 20 | 216 | 0 |

Read-out:
* **k ≤ 2, non-adjacent placements: never an alternative and never a continuum**, for every shape up to
  n = 12 (supports Theorem 2(b) beyond the hand proof).
* **k = 3: alternative ⇔ continuum ⇔ closed claw**, with zero exceptions over all shapes n ≤ 9
  (Theorem 3 + its converse at k = 3).
* **k = 4: alternative ⇔ continuum still holds with zero exceptions**, but closed claws no longer
  explain all cases: the extra ones are *split claws* (two nodes v, w ∈ V(S) at distance 2 through a
  node u ∉ V(S), with the other four neighbours of v and w in V(S); the move α = aδ_v + bδ_w is chosen
  to cancel at u). Since (B) ⇔ (C) is immediate from Theorem 1, the open statement is (A) ⇔ (B).
* Results are identical for unit and random branch lengths: these are combinatorial, not
  "special-length", phenomena.

Explicit generic counterexample (`results/counterexample_k3.txt`, n = 5, lengths 0.7, 0.4, 0.9, 0.5,
0.3, 0.8, 0.6): truth = edges {(2,5),(6,0),(7,1)} with p = (0.2, 0.4, 0.4); the exact LP search finds
{(5,6),(6,0),(7,1)} with p = (0.43, 0.27, 0.31) and {(6,0),(5,7),(7,1)} with p = (0.33, 0.42, 0.25), all
reproducing d exactly. Abundances differ by up to 0.23, not a negligible ambiguity.

### 4b. Parameter level (random parameters; LP over all edge sets)

`code/run_enum.py` + `code/analyze_enum.py`: all labeled topologies for n = 4–6 (n = 7, 8 sampled), every
(or sampled) true edge set with k = 1–3, five parameter regimes (generic random lengths/x/p/ȳ, the same
with ȳ = 0, all-unit-lengths symmetric with x = 1/2 and equal p, symmetric with ȳ = 0, and unit lengths
with random rational p, x). For each noise-free d we search *all* edge sets with |S'| ≤ k+1 and check for
an exact, strictly valid solution (p > 0, 0 < x < 1, ȳ ≥ 0) by least-squares residual + LP. Full table:
`results/enum_summary.md`.

Instances: 1,279 (k=1), 4,495 (k=2), 5,125 (k=3) per regime; 5 regimes ≈ 54k LP-certified instances
(table: `results/enum_summary.md`). Selected rows (generic regime; the other regimes agree except as noted):

| k | true placements | N | continuum on true edges | alternative with k' < k | k' = k | k' = k+1 |
|---|---|---|---|---|---|---|
| 1 | — | 1279 | 0 | 0 | 0 | 0 |
| 2 | non-adjacent | 2997 | 0 | 0 | 0 | 1506 |
| 2 | adjacent | 1498 | **1498 (100%)** | 0 | 0 | 1498 |
| 3 | non-adjacent, no claw | 598 | 0 | 0 | 0 | 583 |
| 3 | non-adjacent, claw (open or closed) | 591 | 530 | 0 | **167 (28%)** | 591 |
| 3 | adjacent | 3936 | 3936 | 1569 | 2749 | 3936 |

* Adjacent placements (allowed by the paper's Assumption 1) always give a continuum, for every
  parameter value (C1). An adjacent triple can even be explained by **fewer** queries (k' = 2 in 1,569
  instances).
* Exact fits with **k+1** placements exist in about half the k = 2 cases and almost all k = 3 cases (C2 and
  "add the edge joining two placements one edge apart"); with ȳ > 0 also k+2 (star split, seen at n = 4–5).
  So "the residual reaches ~0" identifies only the minimal k.
* **Symmetry is not the culprit, contrary to the paper's guess:** with unit lengths, x = 1/2 and equal p,
  claws give a continuum but **0** alternative edge sets (whether a different edge set appears depends
  on which neighbour of the claw centre runs out of mass first, and symmetric values tie). Generic
  values give alternatives in 28% of claw instances.
* How common are the bad configurations? For k uniformly random edges on the E1 trees
  (`results/claw_frequency.md`): an adjacent pair has probability 0.30/0.80/1.00 for k = 5/10/20 on
  birds-jarvis (n = 48) and 0.07/0.30/0.79 on hemipteroid (n = 193); closed claws 0.04/0.30/0.93 and
  0.003/0.03/0.17. In the authors' own simulations (their generator bans queries whose parent or
  grandparent is taken) 63/446 true sets still contain an adjacent pair and 16/446 a closed claw. Their
  exact search "found the true placements in all cases" for k ≤ 3 because edges stay identifiable for
  adjacent pairs; only (p, x) do not.
* Algorithmic by-product: since d determines m up to Laplacian moves and spurious spreading lowers ȳ,
  the single LP **max ȳ s.t. F m + ȳ 1 = d, 1ᵀm = 1, m ≥ 0** (no search, no k) recovers the true node
  measure exactly in 121/121 (k=1), 118/118 (k=2), 98/106, 91/108, 73/98 (k = 3, 4, 5 without claws) random
  instances on trees with n = 6–40 (`results/maxy_lp_noisefree.txt`); it always fails on claws (as it
  must). On noisy data, though, a convex relaxation does not beat DecoDiPhy (§5b).

## 5. Choosing k on the paper's simulation protocol (paired)

Setup (`code/ksel_run.py`): the authors' own generator (`MISC/simulate_exp.py`: random query leaves with
their spacing rule, p uniform-normalized with min p ≥ 0.1/k², noise = moving 10^5 read-pieces an
exponential number of edges, exp_scale 1 ≈ "noise1", 2 ≈ "noise2") on 9 of their 12 E1 biological trees
(n = 30–193), k ∈ {2,3,5,7,10}, 10 seeds, 3 noise levels; then the authors' own greedy search
(`decodiphy.search.hill_climbing`, OSQP small-problem, radius 2) run for **every** k = 1..k_true+3.
Every stopping rule is applied to the *same* trajectory, so comparisons are paired. Free parameters
(BIC penalty λ, elbow ratio τ, min_p) are tuned on 5 trees (bees, birds-jarvis, 1kp, beetles,
hemipteroid) and reported on the other 4 (mammals-song, tilapia, pancrustacean, fish-troyer).
Rules: paper (min_p = 0.01), paper with tuned min_p, BIC n·log(RSS/n) + λ·3k·log n, loss-ratio elbow,
and the theory-motivated **adjacency stop** ("adj": also stop, returning k−1, as soon as the k-solution
contains two placements on adjacent edges — by C1/C2 such pairs are unidentifiable and are what star
splits look like). Jaccard = |chosen edges ∩ true edges| / |union|; Wilcoxon signed-rank on per-instance
Jaccard differences vs the paper's rule.

Runs: 1,350 attempted, **1,338 analysed**. 12 were excluded because the authors' hill climbing never
terminated: it cycles, so a 300 s cap was applied; mammals-song k=7/10 and pancrustacean k=10, all three
noise levels each. Separately, a failure path in their `optimize.py` references an undefined variable
when OSQP fails. Runtime per trajectory (k = 1..k+3): 1–60 s. Full tables: `results/ksel_summary.md`.

**Held-out trees (588 runs; parameters tuned on the other 5 trees):**

| rule | exact-k acc | mean abs(k̂−k) | bias | mean Jaccard | ΔJac vs paper (Wilcoxon p; wins/losses) |
|---|---|---|---|---|---|
| paper (min_p = 0.01) | 0.323 | 1.38 | +0.81 | 0.560 | — |
| paper, min_p tuned (0.005) | 0.338 | 1.57 | +1.38 | 0.566 | +0.006 (0.11; 81/76) |
| BIC (λ tuned = 2) | 0.226 | 2.05 | +2.03 | 0.530 | −0.030 (6e-6; 121/207) |
| loss-ratio elbow (τ tuned = 0.9) | 0.321 | 2.01 | +1.99 | 0.549 | −0.010 (0.07; 116/178) |
| **adjacency stop** (parameter-free) | 0.323 | 1.36 | +0.69 | 0.567 | +0.007 (0.001; 29/10) |
| adjacency stop + tuned min_p | 0.342 | 1.51 | +1.21 | 0.576 | +0.016 (1e-4; 110/78) |
| **learned** (logistic on log min p, relative ȳ drop, log loss ratio, adjacency) | **0.425** | **1.06** | **+0.36** | **0.582** | **+0.022 (2e-6; 131/61)** |
| oracle (true k) | 1 | 0 | 0 | 0.572 | +0.012 (0.89) |

By noise level (held-out trees): noise-free, the paper's rule under-estimates k (true p can be as small
as 0.1/k² < 0.01). Here the loss-ratio rule, i.e. "the first k with ~zero residual", the rule Theorem 1 /
C2 justify, gets 0.964 exact (Jaccard 0.962 vs 0.918). With noise the paper's rule *over*-estimates by about
1.4 (our noise reproduction differs from the paper's Fig. S5, which reports under-estimation). The learned
rule improves noise1 (exact-k 0.153 → 0.270, Jaccard +0.035, p = 7e-4) and is neutral at noise2.

**Leave-one-tree-out over all 9 trees (1,338 runs), learned rule vs paper:** exact-k 0.326 → 0.428
(sign test on runs where exactly one rule is right: 196 vs 60, p = 5e-18), mean abs(k̂−k) 1.38 → 1.11, mean
Jaccard 0.590 → 0.599 (+0.009, Wilcoxon p = 0.003). The gain is driven by noise1 (+0.042, p = 3e-7);
**noise2 is slightly worse (−0.016, p = 0.01)**.

Diagnostics behind the rule, AUC for "this added placement is spurious": min p 0.91/0.92/0.86
(noise 0/1/2); ȳ barely changes for spurious additions while true additions cut it (AUC of −Δȳ
0.93/0.89/0.83); adjacency to an existing placement 0.87/0.63/0.57. Beyond the true k, 42% of newly
added placements (4,412, noisy runs) are adjacent to an existing placement, vs 16% expected for a
random edge: these are C2's star splits.

**Honest caveat: choosing k is not the accuracy bottleneck.** Even the oracle k improves Jaccard by only
+0.018 over the paper's rule on all runs, and is *worse* at noise2 (−0.042): with heavy noise the greedy
search places the extra components wrongly anyway. The k-selection gains are real and significant
but small in placement accuracy (+1–2 Jaccard points).

### 5b. A convex relaxation suggested by Theorem 1 (negative on noisy data)

`code/noisy_qp.py`: fit the node measure directly. Stage 1: non-negative least squares over all
2n−2 nodes. Stage 2: maximize ȳ subject to RSS ≤ (1+τ)·RSS₀. Both are convex (CVXPY/Clarabel); there is
no search over edge sets and no k. Evaluated on the same simulated d̂ as DecoDiPhy, seed 1 (135 paired
runs, all 9 trees × 5 k × 3 noise levels; τ tuned on the training trees but irrelevant in practice). The
metric is the tree earth-mover distance between true and estimated node measures, i.e. weighted
UniFrac on placements with pendant lengths ignored, as in the paper (`results/convex_vs_decodiphy.md`):

| noise | n | DecoDiPhy (paper k) | DecoDiPhy (oracle k) | convex | Δ (convex − DecoDiPhy) | Wilcoxon p | convex better/worse |
|---|---|---|---|---|---|---|---|
| 0 | 45 | 0.0003 | 0.0000 | 0.0005 | +0.0001 | 0.004 | 7/38 |
| 1 | 45 | 0.0185 | 0.0180 | 0.0193 | +0.0008 (+4%) | 4e-5 | 7/38 |
| 2 | 45 | 0.0360 | 0.0379 | 0.0388 | +0.0028 (+8%) | 3e-7 | 8/37 |

Runtime 1.7 s vs 2.8 s per instance. The convex fit spreads mass over 12–15 nodes where the truth has
2k = 10. **Negative result:** the relaxation is not more accurate than DecoDiPhy on noisy data, although
the noise-free LP version is exact for k ≤ 2 (§4b). Making it competitive would need a sparsity
mechanism (reweighting, or projection onto claw-free supports): a possible extension, not a plan.

## 6. Verdict

**Verdict: promising** for a 4-week CS581 project, on the strength of the theory. The k-selection
part is a *measurable but modest* improvement and should be pitched as secondary.

What the pilot established (numbers first):
* A complete algebraic characterization (Theorem 1, fully proved): d determines exactly the node
  measure modulo weighted-Laplacian moves at internal nodes. Proof uses only tree additivity: the
  distance is linear along edges, plus the identity L D + 2I = (2·1 − deg)1ᵀ.
* Exhaustive, parameter-free verification over **all unlabeled shapes up to n = 12** (k ≤ 2) and n ≤ 9
  (k = 3), n ≤ 8 (k = 4); plus about 54k LP-certified random-parameter instances (n = 4–8, k ≤ 3, 5
  regimes). Zero exceptions to: k ≤ 2 non-adjacent ⇒ identifiable; adjacency ⇒ continuum (100%, not a
  corner case); k = 3 alternative ⇔ continuum ⇔ closed claw; alternative ⇔ continuum through k = 4.
* The authors' conjecture that non-identifiability needs "extreme symmetry" is **false**: generic
  lengths give alternatives in 28% of claw instances at k = 3, symmetric parameters in 0%. A concrete
  n = 5 counterexample has abundances differing by 0.23.
* Choosing k on the authors' own simulation protocol (1,338 paired runs, 9 trees): BIC and elbow rules
  are worse than the paper's rule. A theory-motivated adjacency stop is parameter-free, small, and
  significant (+0.007 Jaccard, p = 0.001). A learned 4-feature rule raises exact-k from 0.33 to 0.43
  (leave-one-tree-out, p = 5e-18) and halves the bias, but Jaccard gains are only +0.01–0.02, and even
  the oracle k gives only +0.018. **k is not the main accuracy bottleneck.**
* Negative: a convex node-measure relaxation (Theorem 1) is exact for noise-free k ≤ 2 by LP, but
  slightly *worse* than DecoDiPhy on noisy inputs.

**Weeks 1–4.**
1. Close the proofs: Theorem 2(b) (last case: two T_U-leaves with opposite signs) and the k = 3
   converse (closed claw is necessary), ideally (A) ⇔ (B) in general. Obtain and read the paper's
   supplement SB.3 / Fig. 1d (bioRxiv was behind a captcha from this container), to make sure the
   overlap is only Claims 1–3.
2. Consequences: which functionals of (p, x) are always identifiable? E.g. the mass on any split not
   separating a claw. A noise-tolerance bound for d̂ via the grounded Laplacian L[I, I], analogous to the
   near-additivity results the paper cites. Short algorithmic corollaries: the max-ȳ LP for noise-free
   inputs, and an "identifiability flag" (rank test) on every output.
3. Stopping rule: package the adjacency/learned rule. Re-validate on a second noise model so it is
   not tuned to the authors' simulator (E2-style reads + krepp distances if feasible on 4 cores;
   otherwise a Gaussian-on-d̂ model). Report k accuracy, Jaccard and wUniFrac.
4. Write-up: theorem + proofs + exhaustive tables + k-selection experiment.

**Main risks.** (i) The authors' unread supplement may already contain part of the adjacency analysis
(their code's `check_identifiability` shows awareness). Theorem 1 and the claw results still look new.
(ii) The general (A) ⇔ (B) statement may resist a short proof; the fallback is k ≤ 3 plus exhaustive
evidence. (iii) The k-selection gains are small and tied to the authors' synthetic noise model. Do not
oversell them as an accuracy win: the deliverable is a theory paper with a modest practical corollary.
(iv) Hill-climbing non-termination in DecoDiPhy (12/1350 runs) must be handled with timeouts.
