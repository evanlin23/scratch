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
other, and their simulator never draws such placements — so they are aware adjacency matters, but the
paper treats it as a corner case.)

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
* **C3 (k = 1 is identifiable)** against any alternative with k' ≤ 1 (proof: a T_U-leaf/sibling-pair
  argument, below).

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
have not closed by hand; it is covered by the exhaustive check below (no counterexample, all shapes n ≤ 9).

**Theorem 3 (closed claw, k ≥ 3) — "only if" proved, "if" exhaustive.** Call v ∈ I a *closed claw* of S
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
continuum; (C) ker L[V ∖ V(S), I] ≠ {0} — a one-line rank test. For k ≤ 3, (C) ⇔ closed claw. At k = 4
other patterns appear (e.g. a "double claw" around an internal edge), but (A) ⇔ (B) ⇔ (C) still holds in
every case we enumerated.

## 4. Exhaustive small-case evidence

### 4a. Configuration level (parameter-free; all unlabeled shapes; all matchings S, all S' with |S'| ≤ |S|)

`code/config_enum.py`: by Theorem 1, "S' is an exact alternative for some (p,x,ȳ)" is the feasibility of
a sign pattern on g = −Lα (zero off V(S)∪V(S'), > 0 on V(S')∖V(S), < 0 on V(S)∖V(S')), an LP; checked
for unit and random branch lengths. Full table: `results/config_enum.md`.

RESULTS_CONFIG

### 4b. Parameter level (random parameters; LP over all edge sets)

`code/run_enum.py` + `code/analyze_enum.py`: all labeled topologies for n = 4–6 (n = 7, 8 sampled), every
(or sampled) true edge set with k = 1–3, five parameter regimes (generic random lengths/x/p/ȳ, the same
with ȳ = 0, all-unit-lengths symmetric with x = 1/2 and equal p, symmetric with ȳ = 0, and unit lengths
with random rational p, x). For each noise-free d we search *all* edge sets with |S'| ≤ k+1 and check for
an exact, strictly valid solution (p > 0, 0 < x < 1, ȳ ≥ 0) by least-squares residual + LP. Full table:
`results/enum_summary.md`.

RESULTS_ENUM

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

RESULTS_KSEL

## 6. Verdict

VERDICT
