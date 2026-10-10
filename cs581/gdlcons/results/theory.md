# Theory note: when does ASTRAL-Pro stay consistent under rooting/tagging error? (pure GDL)

*Pilot, 2026-10-10. Branch `claude/cs581-gdlcons`. Every claim below is checked numerically: `code/predict_margin.py` gives exact values and `code/block_se.py` gives simulation estimates with standard errors.*

## 1. Setting and definitions

**Model.**
- Pure GDL with **branch-specific** duplication/loss rates (λ_e, μ_e), as in the sibling session's simulator.
- No ILS, so the gene-family tree is the locus tree, with true root and true D/S event labels.

**ASTRAL-Pro (Zhang et al. 2020, doi:10.1093/molbev/msaa139), from `results/prior_art.md`.**
- *Tagging rule (Def. 1).* A node may be tagged S only if its two children's species sets are disjoint.
- *SQ.* A speciation-driven quartet has four leaves from four species, and the LCA of every three of them is tagged S.
- *Classes.* SQs on the same four species are equivalent iff they share the "anchor LCA". The per-family score of a topology is the number of classes with that topology.
- *Thm 2.* Consistent under GDL "for correctly rooted and tagged gene trees".

**Observation 0.** The unrooted topology of four leaves is fixed by the unrooted gene tree. Rooting and tagging only decide **which** quartets are counted, so errors act as a *selection filter*. With true tags every SQ is orthologous and displays the species quartet (Prop. 3), so the error-free score has no wrong-topology mass at all. Consistency fails only if errors select wrong-topology quartets more often than correct ones.

**Observation 1** (also in the prior-art report). On the **true root** under pure GDL, the species-overlap rule never mislabels a true speciation, because its children have disjoint species sets. Its only errors are *hidden paralogs*: duplications whose two sides end up species-disjoint after complementary losses, which get labelled S. So the realistic tagging error is one-sided (D→S) and only hits hidden paralogs.

## 2. Error models used

| name | rooting | tags | inside ASTRAL-Pro's Def. 1? |
|---|---|---|---|
| `true` | true | true | yes |
| `ovl` | true | species-overlap rule (what ASTRAL-Pro would assign, given the true root) | yes |
| `d2sv(q)` | true | each hidden paralog labelled S with prob. q (q = 1 is `ovl`) | yes |
| `flipv(q)` | true | `d2sv(q)` plus each true S labelled D with prob. q | yes |
| `rovl(p)` | with prob. p, re-root at a uniformly random edge | overlap rule | yes |
| `rsp-X` | root on the leaf edge of a random copy of species X | overlap rule | yes |
| `d2s(q)`, `flip(q)` | true | *any* D labelled S (or any tag flipped) with prob. q, including D nodes whose children share species | **no**. These are S nodes ASTRAL-Pro's own tagger could never produce |
| `own` | ASTRAL-Pro3's own rooting | ASTRAL-Pro3's own tagging (unmodified binary, unrooted input) | yes |

**How rooting and tags are fed in.**
- The patched binary `astral-pro3-fixed` keeps the input root and reads the tags from branch lengths when `APRO_FIXED=1`.
- Otherwise it runs ASTRAL-Pro3's own algorithm unchanged; the patch is about 25 lines in `annotateTree`, documented in `code/README.md`.

## 3. Exact limiting scores for a 4-taxon caterpillar

**Setup.**
- Species tree (((A,B)x,C)y,D), no branch above the root.
- Duplications may occur anywhere. Only those on the **y-branch** (root → y, length T, rates λ, μ) can produce wrong quartets: a wrong quartet needs a duplication above a three-taxon clade.
- D only needs to be present, and its copy number does not matter. Every LCA of a triple that includes d is the gene-tree root, which is S.

**Which nodes can carry a class.** For a quartet (a, b, c, d), the LCA of every triple is either w := LCA(a, b, c) or the root. So a class is counted iff w is tagged S. With the true root:

- **Orthologous classes.** w is the image of the speciation y. That y-copy needs A and B on its x-side and C on its C-side. Expected number per family:
  O = m(T) · P_x(A∧B) · s_C, with m(t) = e^{(λ−μ)t}.
  - Here s_C is the survival probability of one copy entering the C branch.
  - P_x(A∧B) is the probability that a copy entering the x branch leaves descendants in both A and B.
- **Hidden-paralog classes.** w is a duplication on the y-branch at time t whose two subtrees are species-disjoint and together cover {A, B, C}. There are three patterns, and each fixes one topology:
  - ({A,B},{C}) gives AB|CD (correct);
  - ({A,C},{B}) gives AC|BD;
  - ({B,C},{A}) gives AD|BC.
  The two subtrees are i.i.d. copies started at time t. Hence
  H_τ = ∫₀ᵀ 2λ m(t) p_{T−t}(S₁) p_{T−t}(S₂) dt,
  where p_s(S) is the probability that one copy at distance s above y leaves exactly the species set S.
  - p_s is computed exactly from the linear birth–death generating function φ_e. If F_e(T) = P(surviving set ⊆ T), then F_e = φ_e ∘ F_v and F_v = Π_children F_c. Möbius inversion turns these into exact-set probabilities (`code/exact_assoc.py`).
  - Each such node is its own class, and these classes are disjoint from the orthologous ones.

**Limiting per-family ASTRAL-Pro scores.**

| tags | AB\|CD (true) | AC\|BD | AD\|BC |
|---|---|---|---|
| true | O | 0 | 0 |
| `ovl` | O + H_AB | H_AC | H_BC |
| `d2sv(q)` | O + q·H_AB | q·H_AC | q·H_BC |

Filtering families on "≥ 4 species present" divides every entry by the same constant, so the comparison is unaffected. For four taxa, ASTRAL-Pro (exact search) is consistent iff the true column is the strict maximum.

**Proposition (4-taxon caterpillar, pure GDL, correct root).**
1. With true tags, ASTRAL-Pro is consistent whenever O > 0. This is Zhang et al.'s theorem.
2. With species-overlap tags, which is ASTRAL-Pro's own rule given a correct root, it is consistent **iff O + H_AB > max(H_AC, H_BC)**.
3. Under random hidden-paralog mislabelling `d2sv(q)`, it is consistent iff
   q < q* := O / (max(H_AC, H_BC) − H_AB)
   when the denominator is positive, and for all q otherwise.
4. A balanced quartet ((A,B),(C,D)) with no branch above the root has no wrong SQs, so it is always consistent.

A root branch with duplications adds four-species patterns of the same form. I did not work that case out.

**Is H_AC > H_AB possible? Yes.** On the plain single-copy picture with independent survivals, the two pattern probabilities are exactly equal: p(AB)p(C) = p(AC)p(B) = Π_X s_X(1 − s_X). Asymmetry comes from **random, unequal copy numbers on the two sides of a duplication**, created by further duplications on the y-branch, combined with **unequal per-copy survival** in A, B and C. A Monte-Carlo and exact scan of the per-node inequality I1: p(AB)p(C) ≥ max(p(AC)p(B), p(BC)p(A)):
- it fails in 3,580 of 77,606 survivable random rate configurations;
- the containment version I2 fails in 14,423 of 200,000.
So the natural "positive association" proof route does not work in general.

## 4. Counterexample (cand2) and checks

**Configuration.** (((A:1,B:1)x:4,C:4)y:0.5,D:1) with the following rates:

| branch | λ | μ | notes |
|---|---|---|---|
| y | 8 | 0.5 | supercritical, about 42 copies reach y |
| x | 0 | 0 | no events |
| A | 1 | 4 | |
| B | 4 | 8 | |
| C | 0.5 | 1 | |
| D | 0 | 0 | single copy |

Per-copy survival is 3.8% in A, 0.9% in B and 7.3% in C.

**Prediction vs simulation** (`code/predict_margin.py`; values per accepted family, prediction divided by the acceptance rate 0.2265):

| | O | ovl: correct | ovl: AC\|BD | ovl: AD\|BC |
|---|---|---|---|---|
| predicted | 0.00476 | 0.2008 | 0.2563 | 0.2206 |
| direct count of node patterns in 5,000 simulated families | 0.0050 | 0.1944 + O | 0.2628 | 0.2186 |
| ASTRAL-Pro3 scores, 40,000 families (`block_se`) | 0.0049 ± 0.0003 | — | correct − AC\|BD = **−0.060 ± 0.003** (predicted −0.0555) | correct − AD\|BC = −0.017 ± 0.003 (predicted −0.0198) |

`d2sv(0.5)`: correct − AC|BD observed −0.0244 ± 0.0024, predicted −0.0254.

**Species-tree level** (`curve4.py`; fraction of disjoint replicate datasets on which ASTRAL-Pro3 returns a wrong 4-taxon tree):

| tags | 100 families | 500 | 2,000 | 10,000 |
|---|---|---|---|---|
| true | 0/40 | 0/40 | 0/20 | 0/4 |
| ovl | 35/40 | 39/40 | 20/20 | 4/4 |
| **own (unmodified ASTRAL-Pro3, true unrooted gene trees)** | 38/40 | 40/40 | 20/20 | 4/4 |
| ASTRAL-multi (ASTER `astral4` with a gene→species map) | 26/40 | 24/40 | 4/6 | — |

So the unmodified ASTRAL-Pro3 software, given **true** gene trees under pure GDL, converges to the wrong species tree here. The cause is its own species-overlap tagging (hidden paralogs) and not gene-tree error. ASTRAL-multi sits near a three-way tie: block estimate correct − AC|BD = −0.003 ± 0.002. That is not significant, so I make no claim about ASTRAL-multi.

**How common?** From the exact formula on 20,000 random rate configurations (λ, μ ∈ {0, …, 8}, branch lengths 0.05–4; `code/prevalence.py`):
- `ovl` is inconsistent in **12 of 15,434** informative configurations (0.08%);
- 11 of these 12 have a supercritical y-branch, usually λ_y = 8;
- a finite q* exists in 353 configurations, with q* < 1 in only 12.

So this is a corner of parameter space: very high duplication on one short branch above a three-taxon clade, followed by heavy loss. It is not the regime of standard simulations (the DISCO and FastMulRFS data use λ = μ and modest rates).

## 5. Naive random flips (outside Def. 1)

**Configuration cand1.**
- Tree: (((A:1,B:1)x:0.01,C:0.5)y:1,D:1).
- Rates: y (4, 0); x (0.5, 2); A (0.25, 4); B (0.25, 1); C (0.25, 8).
- Here `ovl`, `own` and `d2sv` are all consistent, and the formula predicts consistency for every q.

**Naive model.** Labelling *any* duplication as S with probability q (`d2s`) breaks ASTRAL-Pro already at q ≈ 0.03–0.05. Block estimates (40,000 families):

| q | correct − AD\|BC |
|---|---|
| 0.02 | +0.025 ± 0.002 |
| 0.05 | −0.020 ± 0.003 |
| 0.1 | −0.10 ± 0.006 |

The species tree is wrong in 4/4 datasets of 10,000 families at q = 0.1 (also 4/4 for `flip(0.1)`), and in 0/4 at q = 0.02.

The cause: an S label on a node whose children share species makes ASTRAL-Pro count copy tuples with multiplicity across the two paralogous subtrees. These counts are large (a 7–90× inflation of the total score) and nearly topology-neutral, so a small bias wins. This is an **artifact of an error model that ASTRAL-Pro's own tagger cannot produce**. It shows that "random tag flips" is the wrong formalisation of the slide's question; the error model has to respect Def. 1.

## 6. What a full answer to slide question 1 needs (4-week scope)

1. **Proof.** Write the Proposition properly: the class decomposition, the integral formula, and the conditioning argument. Add the version with a root branch, and the general n-taxon statement. For n > 4, a quartet-wise wrong majority does not automatically mean ASTRAL-Pro is wrong. Use a 4-taxon induced example inside a larger tree, or prove inconsistency via a dominating quartet argument.
2. **Rooting error.** `rovl(p)` and `rsp-X` stayed consistent in every small configuration tried. A formula like the one above, for "random root + overlap tags", is the obvious next step. On random roots, true speciations can be tagged D, which only removes correct mass, and new hidden-paralog-like S nodes appear.
3. **Sufficient conditions.** Find conditions under which O + H_AB > max H_wrong: for example λ_e ≤ μ_e on the branch above every three-taxon clade, or an upper bound on λ_y·T. In the prevalence scan, only 1 of the 12 failures has a subcritical y-branch, so "no supercritical branch" is a natural conjecture. The sibling session found the same condition for its ortholog distance.
4. **ASTRAL-multi.** Check whether Legried et al.'s theorem assumes uniform rates [M], since cand2 sits on a near-tie for ASTRAL-multi.
