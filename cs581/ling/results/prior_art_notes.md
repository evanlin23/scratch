# Prior art: stochastic models of linguistic character evolution and testing tree-estimation methods under them

Literature scan for the CS581 pilot (scan date: 2026-10-09).

**How to read the verification marks**
- **DOI verified** means the DOI resolved in the Crossref API and the returned title, authors, year, volume and pages matched.
- **Read** means I read the full text (PDF converted to text).
- **Abstract/secondary** means I relied on the abstract, a search snippet or my own background knowledge. These are flagged where they occur.

---

## 0. Headline finding

Two directly relevant, openly available resources already exist:

1. **The WERN-2006-style simulator.** Canby et al. 2024 (item 8a) released a Java simulator, **LingPhyloSimulator**. It simulates multi-state characters down trees and networks under two models:
   - the WERN 2006 model (homoplasy, borrowing, deviation from the lexical clock, heterotachy);
   - a new lexical-polymorphism model.

   The repo also includes 32 model trees (30 leaves each), networks with 1–3 contact edges, configs, the simulated datasets, and PAUP*/MrBayes/TraitLab commands.
2. **A corrected and extended Ringe–Taylor IE dataset** (`ie_dataset.csv`, 24 languages, 370 characters, polymorphisms included).

Together these give a 4-week project a ready-made baseline. You can replicate their protocol, then perturb one thing.

---

## 1. Barbançon, Evans, Nakhleh, Ringe & Warnow (2013)

**Citation:** F. Barbançon, S.N. Evans, L. Nakhleh, D. Ringe, T. Warnow. "An experimental study comparing linguistic phylogenetic reconstruction methods." *Diachronica* 30(2):143–170, 2013.
- **DOI:** 10.1075/dia.30.2.01bar (verified)
- **Status:** Read
- **PDF:** http://tandy.cs.illinois.edu/Diachronica-barbanson.pdf; appendix at http://tandy.cs.illinois.edu/diachronica-appendix.pdf
- The Berkeley tech report (732) is an earlier version; a preliminary version appeared at the UCSB "Languages and Genes" meeting.

**Summary.** This is the first large simulation study of linguistic tree estimation under a linguistically motivated model, the WERN 2006 model (item 7). Characters evolve down a genetic tree plus contact edges.

**(a) Simulation design**
- **Data volume:** 3,584 datasets in total. The study ran 28 experimental conditions; each used 32 random model networks with 30 leaves and 0–3 contact edges, with 4 replicate datasets per network.
- **Characters:** 301–360 per dataset. 300 are "lexical" (split into slow, medium and fast); the rest are morphological.
- **Factors varied:** rate of evolution, homoplasy level, deviation from the lexical clock (dlc), heterotachy/deviation from rates-across-sites (het), and number of contact edges, each at low/medium/high.
- **Per-character parameters:** homoplasy factor(c), borrowing(c), dlc(c), het(c).
- **Borrowing:** lexical characters only.
- **Calibration to IE (Nakhleh et al. 2005):**
  - "Screened" setting: 1% of lexical characters homoplastic and 6% borrowed; 0% of morphological characters homoplastic.
  - "Unscreened" setting: 13% of lexical and 24% of morphological characters homoplastic; 7% of lexical characters borrowed.
- **No polymorphism. No missing data.**
- **Methods compared:** UPGMA, NJ, MP, weighted MP, weighted maximum compatibility, and Gray & Atkinson (binary presence/absence encoding analysed in MrBayes under the restriction-site model + gamma).
  - Weighting: morphological and phonological characters got weight 50 as a "proxy for infinity".
  - Accuracy measure: FN and FP bipartition rates against the genetic tree.
- **Findings:**
  - Accuracy ranking: UPGMA < NJ < G&A < MP.
  - Weighting high-confidence (low-homoplasy) characters greatly improves MP and MC on screened data. WMP/WMC should not be used on unscreened data.
  - Screening helps MP more than it helps binary-encoding methods.
  - Every method except UPGMA recovers about 90% of true edges. The differences sit in the "fine details" (e.g., Germanic, Italo-Celtic, Greco-Armenian).

**(b) What it leaves untested (the authors list most of these in §5)**
- Stochastic Dollo and other newer Bayesian methods (Nicholls & Gray 2008; Ryder & Nicholls 2011).
- Network-estimation methods.
- Polymorphic characters. The authors call methods for polymorphism "one of the outstanding problems".
- Longer MCMC runs for G&A; the G&A results may partly be under-converged.
- Separating the effect of binary encoding (multi-state to binary model misspecification) from the effect of the method.
- Missing data.
- ML (rather than Bayesian) on binary or multi-state encodings.
- Ascertainment bias.
- Only 30 taxa, and only 0–3 contact edges.

## 2. Ringe, Warnow & Taylor (2002)

**Citation:** D. Ringe, T. Warnow, A. Taylor. "Indo-European and computational cladistics." *Transactions of the Philological Society* 100(1):59–129, 2002.
- **DOI:** 10.1111/1467-968X.00091 (verified)
- **PDF:** http://www.cs.rice.edu/~nakhleh/CPHL/RWT02.pdf (HTTP 200)

**Summary.**
- Introduces the curated 24-language IE dataset (lexical, morphological and phonological characters, multi-state).
- Uses a weighted maximum compatibility / perfect-phylogeny approach. Hand-screening for homoplasy and borrowing gives high weight to morphological and phonological characters.
- Proposes the RWT IE tree. Germanic appears non-treelike, which motivated the later network work.

**(a) Simulation:** None. It is an empirical analysis.

**(b) What it leaves untested**
- Method accuracy (the true tree is unknown).
- Formal statistical modelling.
- Handling of polymorphic characters, which were removed or split-coded.

## 3a. Nakhleh, Ringe & Warnow (2005), *Language*

**Citation:** L. Nakhleh, D. Ringe, T. Warnow. "Perfect phylogenetic networks: a new methodology for reconstructing the evolutionary history of natural languages." *Language* 81(2):382–420, 2005.
- **DOI:** 10.1353/lan.2005.0078 (verified)
- **Status:** Partially read. PDF: http://www.cs.rice.edu/~nakhleh/Papers/81.2nakhleh.pdf

**Summary.**
- Extends perfect phylogeny to a genetic tree plus bidirectional "contact edges". A character is compatible with the network if it is compatible with at least one tree displayed by the network.
- Formulates the problem as Minimum Increment to Perfect Phylogenetic Network (MIPPN), solved by the MIPPN software (listed on the CPHL page, but not visibly downloadable).
- On the screened IE data, a few contact edges (mainly involving Germanic) make almost all characters compatible.

**(a) Simulation:** None. Empirical IE only.

**(b) What it leaves untested**
- Accuracy of MIPPN at detecting true borrowing on simulated data.
- Whether contact edges can be distinguished from homoplasy.
- Behaviour on polymorphic data.

## 3b. Nakhleh, Warnow, Ringe & Evans (2005), *TPS*

**Citation:** L. Nakhleh, T. Warnow, D. Ringe, S.N. Evans. "A comparison of phylogenetic reconstruction methods on an Indo-European dataset." *Transactions of the Philological Society* 103(2):171–192, 2005.
- **DOI:** 10.1111/j.1467-968X.2005.00149.x (verified)
- **Status:** Read the opening. PDF: http://www.cs.rice.edu/~nakhleh/Papers/TPS05.pdf (Berkeley TR 672 has the full version).

**Summary.**
- Runs UPGMA, NJ, MP, weighted MP, MC/WMC and G&A on 24 IE languages × 336 characters (screened and unscreened; full vs lexical-only).
- The methods agree on the major subgroups and differ on Germanic and other fine structure.
- Lexical-only analyses look less reliable. Some relationships appear only when morphological and phonological characters are included.

**(a) Simulation:** None. This empirical study motivated the 2013 simulation paper.

**(b) What it leaves untested:** Which method is actually more accurate, since there is no ground truth.

## 4. Bouckaert et al. (2012), *Science*

**Citation:** R. Bouckaert, P. Lemey, M. Dunn, S.J. Greenhill, A.V. Alekseyenko, A.J. Drummond, R.D. Gray, M.A. Suchard, Q.D. Atkinson. "Mapping the origins and expansion of the Indo-European language family." *Science* 337(6097):957–960, 2012.
- **DOI:** 10.1126/science.1219669 (verified)
- **Status:** Abstract and background knowledge only; not re-read.

**Summary (from memory, partially unverified).**
- Bayesian phylogeography in BEAST: binary cognate presence/absence, CTMC/covarion and stochastic Dollo variants, relaxed clock, and a relaxed random-walk spatial diffusion model.
- Data: about 103 ancient and modern IE languages × about 200 meanings (Swadesh-type, IELex-derived).
- Supports an Anatolian origin with a root around 8,000–9,500 BP.

**(a) Simulation:** No tree-accuracy simulation study that I know of. The supplement contains model comparisons and robustness checks (unverified).

**(b) What it leaves untested**
- Accuracy under borrowing and polymorphism.
- Effects of binary encoding.
- The impact of ancestry constraints; Chang et al. 2015 later argued these matter.

## 5. Chang, Cathcart, Hall & Garrett (2015), *Language*

**Citation:** W. Chang, C. Cathcart, D. Hall, A. Garrett. "Ancestry-constrained phylogenetic analysis supports the Indo-European steppe hypothesis." *Language* 91(1):194–244, 2015.
- **DOI:** 10.1353/lan.2015.0005 (verified)
- **Status:** Background knowledge only; not re-read.

**Summary.**
- Re-analyses IELex-derived binary cognate data in BEAST.
- Constrains attested ancient languages (e.g., Latin, Old Irish, Vedic) to be direct ancestors of their modern descendants. Without the constraints, branch lengths come out too long.
- Discusses "advergence": independent parallel semantic shifts, i.e. homoplasy from shared semantic change, which inflates the age estimates.
- The constrained analyses give a root of about 6,500 BP (steppe).

**(a) Simulation:** I do not recall a full tree-accuracy simulation (unverified).

**(b) What it leaves untested:** Topology accuracy as a function of homoplasy, polymorphism or borrowing; the focus is dating.

## 6. Stochastic Dollo family (Nicholls, Gray, Ryder, Kelly)

**6a.** G.K. Nicholls, R.D. Gray. "Dated ancestral trees from binary trait data and their application to the diversification of languages." *Journal of the Royal Statistical Society: Series B* 70(3):545–566, 2008.
- **DOI:** 10.1111/j.1467-9868.2007.00648.x (verified)
- **Model:** Stochastic Dollo (SD) for binary cognate traits. Each trait is born once (no homoplasy) and dies at rate μ. The model adds catastrophes (rate-heterogeneity bursts) and an ascertainment correction.
- **Data:** Applied to IE. I recall that this includes the Ringe–Taylor data plus Dyen et al.; this is unverified.
- **Validation:** I recall checks on synthetic data (unverified).
- **Untested:** Robustness to borrowing and homoplasy. Multi-state structure is ignored, since cognate classes within a meaning are treated as independent.

**6b.** R.J. Ryder, G.K. Nicholls. "Missing data in a stochastic Dollo model for binary trait data, and its application to the dating of Proto-Indo-European." *Journal of the Royal Statistical Society: Series C* 60(1):71–92, 2011 (online 2010).
- **DOI:** 10.1111/j.1467-9876.2010.00743.x (verified)
- **Contribution:** Integrates over missing-data patterns in SD and fits IE.
- **Validation:** Model checking on synthetic data (unverified in detail).
- **Untested:** Borrowing; polymorphism as a process.

**6c.** L.J. Kelly, G.K. Nicholls. "Lateral transfer in Stochastic Dollo models." *Annals of Applied Statistics* 11(2):1146–1168, 2017.
- **DOI:** 10.1214/17-AOAS1040 (verified; the DOI in the brief was wrong, since 10.1214/17-AOAS1028 is an unrelated paper). arXiv: 1601.07931.
- **Model:** SD plus lateral transfer. Traits are copied between contemporaneous lineages at a rate. The likelihood requires solving large ODE systems.
- **Data:** Applied to Eastern Polynesian lexical data, where it fits better than SD without borrowing.
- **Simulation:** Synthetic-data validation of parameter recovery (abstract and secondary sources; detail unverified).
- **Untested:**
  - A head-to-head topology-accuracy comparison against parsimony or compatibility under a non-Dollo generating model (e.g., WERN).
  - Scalability beyond small families.

**TraitLab:** a Matlab package for SD with catastrophes, missing data and lateral transfer, including a simulator. It is described in arXiv 2308.09060 (2023). The old URL https://www.stats.ox.ac.uk/~nicholls/TraitLab/ returned 404; get it via the arXiv paper or GitHub.

## 7. Warnow, Evans, Ringe & Nakhleh (2006): the WERN model

**Citation:** T. Warnow, S.N. Evans, D. Ringe, L. Nakhleh. "A stochastic model of language evolution that incorporates homoplasy and borrowing." In P. Forster & C. Renfrew (eds.), *Phylogenetic Methods and the Prehistory of Languages*, pp. 75–87. McDonald Institute for Archaeological Research, Cambridge, 2006.
- **DOI:** None found; book chapters from this series often lack one. The page numbers come from memory and are **unverified**.
- **Status:** Read the Berkeley TR 673 version: https://stat.berkeley.edu/sites/default/files/tech-reports/673.pdf

**Summary.**
- Multi-state characters with an unbounded state space. Each character has a small set of "homoplastic states" that can arise more than once; all other changes produce novel states.
- The genetic tree is augmented with contact edges, along which a character can be borrowed.
- Allows character-specific rates: deviation from the lexical clock and heterotachy.
- **Proves identifiability of the tree under no borrowing, provided homoplastic states are known in advance.** It gives statistically consistent methods and linear-time likelihoods under these assumptions.

**(a) Simulation:** None in the chapter. It is the generating model used in Barbançon 2013 and (as `CharacterClass`) in Canby 2024.

**(b) What it leaves untested**
- Polymorphism.
- Identifiability with borrowing. The 2023 follow-up shows the network topology is identifiable under mild constraints (T. Warnow, S.N. Evans, L. Nakhleh, "Progress on constructing phylogenetic networks for languages", arXiv:2306.06298, in *The Method Works: Studies on Language Change in Honor of Don Ringe*, Springer 2023, pp. 45–62; chapter DOI unverified).
- The 2023 network algorithms are theoretical; no implementation or simulation is reported.
- Unknown homoplastic states.

## 8. Post-2009 and 2018–2026 work

### 8a. Canby, Evans, Ringe & Warnow (2024): the most important recent item

**Citation:** M.E. Canby, S.N. Evans, D. Ringe, T. Warnow. "Addressing polymorphism in linguistic phylogenetics." *Transactions of the Philological Society* 122(2):191–222, 2024.
- **DOI:** 10.1111/1467-968X.12289 (verified)
- **Status:** Read. PDF: http://tandy.cs.illinois.edu/Canby-Transactions2024.pdf
- **Code and data:** https://github.com/marccanby/LingPhyloSimulator (public; last commit Feb 2026)

**Summary.**
- New multi-state model in which lexical change happens only through polymorphism: a new state arises through semantic shift or borrowing, coexists with the old one, and then the old one dies. Per-state death rate i·μ limits how large polymorphism can grow.
- Re-analysis of IE with the polymorphic characters restored. Two equally optimal MP trees are found; one matches Nakhleh et al. 2005, the other moves Italo-Celtic relative to Tocharian.

**(a) Simulation design**
- **Data volume:** 20 conditions × 128 datasets = 2,560 datasets. The 32 trees (30 leaves) and the replicate scheme follow Barbançon 2013.
- **Conditions:** 5 polymorphism levels (16/33/50/67/80% of characters allowed to be polymorphic) × 0–3 contact edges.
- **Characters per dataset:** 300 lexical (slow/medium/fast, 13% homoplastic) and 20 morphological (25% homoplastic, not polymorphic). Moderate dlc and heterotachy.
- **Calibration:** Parameters were set from the *unscreened* IE data, which the authors note is deliberately not favourable to MP. About 67% of IE lexical characters are polymorphic, which is the "high" level.
- **No missing data.**
- **Methods compared:**
  - MP variants 1–4 on polymorphic multi-state characters (MP4 is a simple modification of the parsimony criterion).
  - MP on the binary encoding; Dollo parsimony.
  - UPGMA and NJ.
  - GA2003 (MrBayes on binary characters).
  - Stochastic Dollo (TraitLab).
  - Maximum compatibility is described in the text; I did not confirm whether it was benchmarked.
- **Findings:**
  - MP on multi-state polymorphic characters (MP4) is best in almost all conditions.
  - GA2003 is second, and matches MP only with no borrowing and high or very-high polymorphism.
  - SD and the distance methods are clearly worse.
  - Bayesian analyses of binary encodings are "less accurate than MP ... especially – but not only – when there is borrowing".

**(b) What it leaves untested (from the authors' future-work section and my reading)**
- Likelihood or Bayesian inference *under their own polymorphism model*; MP is the only model-aware method.
- Mostly-modern-language settings, where borrowing and polymorphism are higher.
- Newer Bayesian models: covarion, BEAST2 multistate, contacTrees.
- Ascertainment-bias correction after character removal.
- Dating.
- Network estimation; MIPPN is not extended to polymorphism.
- Missing data.
- Weighted compatibility / bounded-homoplasy methods on polymorphic data.
- Borrowing-detection filters.
- ML on multi-state characters (IQ-TREE/RAxML-NG MK/GTR-type multistate models), as opposed to binary encodings.

### 8b. Greenhill, Currie & Gray (2009)

**Citation:** S.J. Greenhill, T.E. Currie, R.D. Gray. "Does horizontal transmission invalidate cultural phylogenies?" *Proceedings of the Royal Society B* 276:2299–2306, 2009.
- **DOI:** 10.1098/rspb.2008.1944 (verified)
- **Status:** PMC full text summarised (PMC2677599).

**(a) Simulation design**
- Binary cognate traits evolved under stochastic Dollo (TraitLab) on 2 fixed "true" trees: one balanced, one unbalanced and chain-like. Root ages are about 5,000 years.
- Borrowing scenarios:
  - local: only between lineages whose common ancestor is within 1,000 or 3,000 years;
  - global: between any contemporaneous lineages.
- Borrowing rates range from 0 to about 50% of traits per 1,000 years.
- Inference used TraitLab SD *without* borrowing; accuracy was measured by quartet distance and root age.

**Findings**
- Topology is robust to "realistic" local borrowing: under 1% perturbation for the balanced tree and about 7.6% for the unbalanced tree.
- Global borrowing degrades the topology roughly linearly, reaching 53–57% at the worst setting.
- Dates are underestimated under borrowing.

**(b) What it leaves untested**
- Only 2 trees, with limited replication.
- The generating and inference models match except for borrowing, so the test favours the method.
- No homoplasy, no polymorphism, no multi-state characters.
- No comparison with parsimony or compatibility.
- No systematic (directional, wholesale) borrowing.

Neureiter et al. 2022 argue that its continuous single-word borrowing is unrealistic.

### 8c. Neureiter et al. (2022): contacTrees

**Citation:** N. Neureiter, P. Ranacher, N. Efrat-Kowalsky, G.A. Kaiping, R. Weibel, P. Widmer, R.R. Bouckaert. "Detecting contact in language trees: a Bayesian phylogenetic model with horizontal transfer." *Humanities and Social Sciences Communications* 9:205, 2022.
- **DOI:** 10.1057/s41599-022-01211-7 (verified)
- **Status:** Abstract and secondary sources only; the full text was blocked.
- **Code:** BEAST2 package contacTrees; IE case study on Zenodo, record 6563028.

**Summary.**
- Models contact as a small number of discrete contact edges, at which each word is borrowed with some probability. This is close in spirit to WERN's contact edges, but in a binary/Bayesian setting.
- Simulations show that it recovers the simulated contact events.
- Ignoring contact biases tree height, rates and topology.
- On IE, it recovers known loans and gives a tree height closer to accepted estimates.

**(b) Untested:** I did not find a comparison against MP or compatibility, or against non-binary generating models. The simulation parameters are unverified.

### 8d. King (2026)

**Citation:** B. King. "Testing the validity and adequacy of linguistic phylogenetic analyses." *PLOS Computational Biology* 22(5):e1014312, 2026.
- **DOI:** 10.1371/journal.pcbi.1014312 (verified via PubMed, PMID 42160374)

**Summary.**
- Simulation-based calibration of the standard BEAST2 phylolinguistic set-up: binary characters with ascertainment correction, meaning-partitions, covarion model, UCLN clock.
- Reweighting partition rates by cognate count (the default) gives poorly calibrated posteriors.
- Posterior predictive checks show the covarion model misfits real lexical data, "likely due to the prevalence of semantic shift and the non-independence of cognate substitutions".

**(b) Untested:** Topology accuracy versus non-Bayesian methods; borrowing.

This is strong independent support for the "binary encoding is misspecified" theme.

### 8e. Other items, briefly

- **Heggarty et al. 2023.** P. Heggarty, C. Anderson, M. Scarborough, B. King, R. Bouckaert, L. Jocz, et al. "Language trees with sampled ancestors support a hybrid model for the origin of Indo-European languages." *Science* 381:eabg0818, 2023. DOI 10.1126/science.abg0818 (verified).
  - Introduces the IE-CoR database (about 160 languages × 170 meanings).
  - Bayesian analysis with sampled ancestors; root about 8,100 BP.
  - No tree-accuracy simulation that I know of (unverified).
- **Kolipakam et al. 2018.** V. Kolipakam, F.M. Jordan, M. Dunn, S.J. Greenhill, R. Bouckaert, R.D. Gray, A. Verkerk. "A Bayesian phylogenetic study of the Dravidian language family." *Royal Society Open Science* 5:171504. DOI 10.1098/rsos.171504 (verified).
  - Empirical BEAST analysis of 20 Dravidian languages; root about 4,500 BP.
  - No simulation.
- **Hoffmann et al. 2021.** K. Hoffmann, R. Bouckaert, S.J. Greenhill, D. Kühnert. "Bayesian phylogenetic analysis of linguistic data using BEAST." *Journal of Language Evolution* 6(2):119–135. DOI 10.1093/jole/lzab005 (verified).
  - Tutorial and review of BEAST2 models for language data.
  - No method-accuracy simulation.
- **Rama & List 2019.** T. Rama, J.-M. List. "An automated framework for fast cognate detection and Bayesian phylogenetic inference in computational historical linguistics." *ACL 2019*, pp. 6225–6235. DOI 10.18653/v1/P19-1627 (verified; P19-1597 is a different paper).
  - LexStat/SCA-style cognate detection (BipSkip) plus MAP tree search (MAPLE).
  - Evaluated against gold cognates and expert trees, not simulation.
  - Related: T. Rama, J.-M. List, J. Wahle, G. Jäger, "Are automatic methods for cognate detection good enough for phylogenetic reconstruction in historical linguistics?", *NAACL 2018*, pp. 393–400, DOI 10.18653/v1/N18-2063 (verified).
- **Jäger 2018.** G. Jäger. "Global-scale phylogenetic linguistic inference from lexical resources." *Scientific Data* 5:180189. DOI 10.1038/sdata.2018.189 (verified).
  - Automated ASJP-based character extraction and ML trees at global scale.
  - Evaluated against Glottolog, not simulation.
- **List et al. 2022.** J.-M. List, R. Forkel, S.J. Greenhill, C. Rzymski, J. Englisch, R.D. Gray. "Lexibank, a public repository of standardized wordlists with computed phonological and lexical features." *Scientific Data* 9:316. DOI 10.1038/s41597-022-01432-0 (verified).
  - Data infrastructure (CLDF), no simulation.
- **Cathcart 2018.** C.A. Cathcart. "Modeling linguistic evolution: a look under the hood." *Linguistics Vanguard* 4(1). DOI 10.1515/lingvan-2017-0043 (verified).
  - Critical discussion of model assumptions in Bayesian phylolinguistics (semantic shift, independence). I did not confirm whether it contains a simulation.
- **Nichols & Warnow 2008.** J. Nichols, T. Warnow. "Tutorial on computational linguistic phylogeny." *Language and Linguistics Compass* 2(5):760–820. DOI 10.1111/j.1749-818X.2008.00082.x (verified).
  - Survey of methods and datasets.
- **Bonet, Phillips, Warnow & Yooseph 1999.** M. Bonet, C. Phillips, T. Warnow, S. Yooseph. "Constructing evolutionary trees in the presence of polymorphic characters." *SIAM Journal on Computing* 29(1):103–131. DOI 10.1137/S0097539796324636 (verified).
  - Combinatorial complexity of polymorphic perfect phylogeny; the model is semantic shift only, with no homoplasy.
- **Splits networks on language data.** Bryant, Filimon & Gray 2005, "Untangling our past: languages, trees, splits and networks", in Mace, Holden & Shennan (eds.), *The Evolution of Cultural Diversity*. Not re-checked, no DOI. NeighborNet is widely used descriptively; I found **no linguistic simulation study benchmarking NeighborNet or other network methods against known borrowing**.
- **Pseudo-Dollo models.** Bouckaert & Robbeets, bioRxiv 10.1101/207571 (2017/18; cited by Canby et al.). A 3-state CTMC approximation to Dollo that tolerates borrowing and coding errors better.

**Searched for but not found (as of 2026-10):**
- A 2022–2026 study comparing *compatibility* methods with Bayesian methods on simulated borrowing data, other than Canby 2024.
- Any "Kim … Warnow" or "Zhang … Warnow" linguistic polymorphism paper. The relevant Warnow-group work is Canby et al. 2024 and Warnow/Evans/Nakhleh 2023.
- Any GQ/quartet-based linguistic method study.

---

## 9. Dataset access (all fetched 2026-10-09)

### CPHL datasets (Ringe & Taylor; 24 IE languages)

The project page **http://tandy.cs.illinois.edu/histling.html** (HTTP 200) lists the datasets. Its Datasets section links to files hosted at Rice; `web.engr.illinois.edu/~warnow/...` redirects to Shibboleth login, and `https://www.cs.utexas.edu/~tandy/histling.html` also returns 200. The Rice index http://www.cs.rice.edu/~nakhleh/CPHL/ is also reachable.

**Screened dataset:** http://www.cs.rice.edu/~nakhleh/CPHL/IEDATA_112603 (200; 19 KB; plain ASCII)
- Format: one header line of 24 language codes (`HI AR GK AL TB VE AV OC LI OE OI LA LU LY TA PE PR LT GO ON OG WE OS UM`).
- Then **294 character lines**, each `id name s1 … s24` with integer states. 22 phonological (P*), 13 morphological (M*), 259 lexical.
- Polymorphic meanings appear as multiple split-coded characters (e.g., `all1`, `all2` with ids `1.1`, `1.2`).
- 30 lines begin with `!`. I believe these mark characters to exclude or treat specially (e.g., uninformative or in some analyses down-weighted), but **I could not confirm the meaning**. Check the coding PDFs or ask.
- There are no explicit missing-data symbols; unknowns appear to be coded as unique singleton states.

**Unscreened dataset:** http://www.cs.rice.edu/~nakhleh/CPHL/IEDATA_050704 (200; 376 lines = header + 375 characters; same format; 30 `!` lines).

**Documentation (all PDFs, all HTTP 200):**
- Phonological coding: http://www.cs.rice.edu/~nakhleh/CPHL/code-p-07.pdf
- Morphological coding: http://www.cs.rice.edu/~nakhleh/CPHL/code-m-07.pdf
- Wordlists: http://www.cs.rice.edu/~nakhleh/CPHL/ie-wordlist-07.pdf
- Lexical coding: http://tandy.cs.illinois.edu/cognations-2k-revised.pdf

**2024 corrected and extended Ringe & Taylor dataset (with polymorphism):** https://raw.githubusercontent.com/marccanby/LingPhyloSimulator/main/example/ie_dataset.csv (200)
- CSV with columns `id,feature,weight,HI,…,UM`; 370 characters.
- The `weight` column is 10000 for phonological characters.
- Polymorphic cells are written `a/b` (224 lines contain one); empty cells mean missing.
- The repo can be cloned with `git clone https://github.com/marccanby/LingPhyloSimulator`, which worked. In this sandbox the github.com web page returned 403 through the proxy, but git worked.
- The repo also has `Simulator.java`, `example/trees.txt` (32 trees, 30 taxa), `example/networks/`, `example/configs/` and `example/simulated_data/` (10 condition folders, e.g., `high_borrowing` with 384 CSVs), plus `software_commands.pdf` (PAUP*/MrBayes/TraitLab commands) and `algorithmic_description.pdf`.

**Old Barbançon simulator:** The CPHL page's links (http://www.cs.utexas.edu/users/francois/logos/Simulator.tar.gz and its documentation PDF) return **404**. LingPhyloSimulator re-implements it (`CharacterClass` = WERN 2006).

### IE-CoR (Heggarty et al. 2023)

- Web: https://iecor.clld.org/ (200). The current release is v1.2.1.
- CLDF dataset on Zenodo: https://zenodo.org/records/23183358 (DOI 10.5281/zenodo.23183358, published 2026-10-06; the older link zenodo.org/records/8089433 redirects there).
- Zip: https://zenodo.org/records/23183358/files/lexibank/iecor-v1.2.1.zip?download=1 (200, about 6.1 MB).
- Data: binary-codable cognate sets per meaning, about 160 languages × 170 meanings. Polymorphism shows up as multiple cognate sets per meaning per language.

### Lexibank and IELex

- Lexibank: https://lexibank.clld.org/ (200) and the Zenodo community https://zenodo.org/communities/lexibank (200). The datasets are CLDF wordlists, many with expert cognates.
- The GitHub organisation pages for lexibank returned 403 through the sandbox proxy; this was not verified further.
- IELex (ielex.mpi.nl) failed with an expired TLS certificate. IELex is effectively superseded by IE-CoR.

---

## What remains untested

These are gaps where a 4-week simulation project could add evidence. All can be built on LingPhyloSimulator and its 32 trees and network configs.

1. **ML on multi-state characters vs. binary encoding vs. MP under the Canby polymorphism model.**
   - Canby et al. tested only MrBayes (binary) and TraitLab SD as model-based methods.
   - Options: IQ-TREE or RAxML-NG with an MK/GTR-type multistate model on the monomorphic or majority-state recoding, and binary + ASC ML.
   - This directly separates "statistical vs. parsimony" from "binary-encoding misspecification".
2. **Treatment of polymorphism as an input transformation.**
   - Compare: drop polymorphic characters (the 2002/2005 practice); split coding; random-state or most-frequent-state resolution; binary encoding; MP4 directly.
   - No study isolates how much accuracy is lost by each transformation, holding the method fixed.
3. **Borrowing-detection filters before tree estimation.**
   - Screen characters by incompatibility count or conflict with a quick tree, then re-estimate.
   - Barbançon 2013 compares only "screened" data generated with fewer homoplastic characters by construction. No one has tested a *data-driven* screening rule on simulated data where false-positive and false-negative screening is possible.
4. **Character-type weighting when weights are mis-specified.**
   - Barbançon showed that WMP/WMC help when high-weight characters really are low-homoplasy, and hurt on unscreened data.
   - The dose-response is untested: what happens as the fraction of homoplastic "high-weight" morphological characters rises from 0 to 25%? And with the 2024 dataset's 10000 weights?
5. **Compatibility with bounded homoplasy, and on polymorphic data.**
   - Maximum compatibility is the classic Ringe–Warnow method, but no simulation evaluates it on polymorphic characters.
   - Also untested: "k-homoplasy" variants (allowing up to k homoplastic events per character) against the WERN model's known homoplastic states.
6. **Stress-testing the model calibration itself.**
   - All the Warnow-group simulations use 30 taxa, 0–3 contact edges, 320–360 characters and parameters tuned to IE.
   - Untested: more contact edges, or clustered/systematic borrowing (Greenhill et al. flag this gap too); more taxa (50–100); fewer characters (Swadesh-100/200 only, the lexical-only setting).
7. **Missing data.** None of the Barbançon or Canby simulations has missing data, yet real IE and IE-CoR data do. Even a simple MCAR or taxon-specific missingness experiment would be new for these models.
8. **Network-aware or borrowing-aware estimators vs. tree estimators on WERN/Canby data.**
   - No study runs contacTrees, TraitLab with lateral transfer, NeighborNet or MIPPN on data generated by a *different* (multi-state, homoplastic) model, and scores both the tree and the recovered contact edges.
   - A small version (few datasets, short MCMC) could be feasible; an FN/FP metric for contact edges is easy to compute.
9. **Model adequacy across generators.**
   - King 2026 shows that covarion misfits real data. Canby shows that binary Bayesian methods lose accuracy on polymorphic data.
   - No one has asked which generator (WERN, Canby-polymorphic, SD + lateral transfer) best reproduces summary statistics of real IE data (e.g., incompatibility distribution, polymorphism rate, consistency index).
   - A posterior-predictive-style comparison against the 2024 Ringe & Taylor and IE-CoR data would support choosing a "realistic" simulation condition.
10. **Heterotachy and lexical-clock deviation × polymorphism interaction.** Canby fixed dlc and het at "moderate". Crossing these with polymorphism and borrowing levels is untested.

## Caveats and honesty notes

- **Read in full or in large part:** Barbançon 2013, Canby 2024, the WERN TR, the 2023 arXiv paper, and the openings of Language 2005 and TPS 2005.
- **Abstracts, PMC summaries or memory only:** Bouckaert 2012, Chang 2015, Nicholls & Gray 2008, Ryder & Nicholls 2011, Kelly & Nicholls 2017, Neureiter 2022, Heggarty 2023, Kolipakam 2018, Hoffmann 2021, Cathcart 2018.
  - Their "(a) simulation" statements are marked unverified where relevant.
- **Page numbers** for the WERN 2006 chapter (75–87) are from memory and unverified.
- **The `!` prefix** in the CPHL data files is uninterpreted; check the coding PDFs or ask.
- **The brief's DOIs for Kelly & Nicholls and for Rama & List were wrong.** The corrected DOIs above are verified.
