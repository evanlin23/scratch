# DL subset trees + GTM: does "deep learning on small subsets, merged by a disjoint tree merger" give accurate 1000-taxon trees?

CS581 (Warnow, Fall 2026) project pilot, about 4.5 h wall-clock on a 4-core, 15 GB, CPU-only machine.

- Code: `cs581/dldtm/code/`. Result tables: `cs581/dldtm/results/`; `SUMMARY.md` has every table with paired statistics.
- Per-replicate trees and logs stayed outside the repository in `/opt/dldtm_runs`. They can be regenerated with the scripts below.

## 0. Verdict

**Not promising** as framed: "DL subset trees + GTM beat the published DTM pipeline or full-data ML."

With the public pretrained models, DL subset trees are much worse than the cheapest classical subset trees. GTM cannot repair that.

**Headline numbers** (FN %, mean over replicates):

| Method | ROSE | RNASim1000 |
|---|---|---|
| GTM + NeuralNJ subsets (≤50 taxa) | 18.1–19.7 | 23.5 |
| GTM + JC-distance NJ subsets | 14.8–16.0 | — |
| GTM + IQ-TREE subsets | 10.1–11.0 | 14.4 |
| FastTree, full data | 9.3–9.8 | 14.2 |
| IQ-TREE, full data | 9.5–10.2 | 15.3 |

**Paired results** (11 replicates pooled):

| Comparison | Mean difference (FN points) | W/T/L | Wilcoxon p |
|---|---|---|---|
| GTM+NeuralNJ vs GTM+IQ-TREE | +8.3 | 0/0/11 | 0.001 |
| GTM+NeuralNJ vs FastTree full | +9.3 | 0/0/11 | 0.001 |

**Phyloformer** (protein weights fed DNA): GTM+PF 16.9–18.8% on ROSE R0, worse than every baseline in 3/3 replicates.

**The cause is domain shift, not the merger and not our harness:**
- On NeuralNJ's own released test data our CPU harness reproduces the paper. NeuralNJ ≈ IQ-TREE (`--fast`) and is better than FastTree (§4.3).
- On ROSE subsets, NeuralNJ stays bad when gap columns are removed and when the subsets are re-simulated under GTR+G on the true subtrees (§4.4). The shift is in the tree and branch-length regime, not only in indels.

**What would be needed:** a "DL + DTM" project is only viable as a *training* project. A DL model would have to be (re)trained or fine-tuned on subtrees that look like centroid-decomposition subsets of large trees, and on CPU that is a risky 4-week job (§6). The more defensible deliverable is a negative or diagnostic study, "DL tree estimators do not transfer to DTM subsets." It would extend Zaharias, Grosshauser & Warnow 2022 from quartets to 50-taxon subsets.

## 1. Question and source

- **Slide source:** first-day lecture (`notes/lectures/CS581-Fall2026-firstday.txt`, "Open problems (and possible course projects)"), item "Deep learning to improve large-scale phylogenetic tree and network estimation."
- **Related slide:** the divide-and-conquer deck (`581-Divide-and-Conquer-trees-2023.txt`) says DTMs "are generic methods, that can be used with any phylogeny estimation method."
- **Pipeline from that deck:** starting tree → centroid decomposition → subset trees → GTM.
- **The idea:** DL tree estimators are limited to roughly ≤100–200 taxa, so let a DTM do the scaling.

## 2. Prior art (details and DOIs in `results/prior_art.md`)

- **Phyloformer**: Nesterenko et al., *MBE* 2025, doi:10.1093/molbev/msaf051; github.com/lucanest/Phyloformer.
  - Public weights (~1.5 MB, 308k parameters), runs on CPU.
  - **Protein only** (alphabet `ARNDCQEGHILKMFPSTWYVX-`).
  - The code hard-caps at 200 sequences (`SEQ2PAIR = seq2pair(200)`). Memory grows as O(n²·L).
  - Outputs a distance matrix; FastME builds the tree.
  - No nucleotide model exists. Phyloformer 2 (arXiv:2510.12976) is also protein-only, and its conclusion lists supertree or ensemble scaling as future work.
- **NeuralNJ** (the "Zhang et al. MBE 2025" reading): doi:10.1093/molbev/msaf260; github.com/ZhangXinru99/NeuralNJ.
  - **DNA** (GTR+I+G).
  - One public checkpoint (5.3 MB), trained on 50-taxon AliSim data with 128–1024 sites.
  - Tested up to about 100 taxa on GPU.
  - It needs a C++ `raxmlpy` binding. That binding failed to build against current raxml-ng (a coraxlib `-fPIC` link error). Plain argmax inference without branch-length optimization needs only the binding's pure-Python helpers, so we use a shim (`/opt/src/raxmlpy_shim`) and a Python 3.10 env, because ete3 breaks on Python 3.13.
- **Fusang**: Wang et al., *NAR* 2023, doi:10.1093/nar/gkad805. DNA, CPU, MIT license, 4–40 taxa.
  - **The only published DL + DTM pipeline we found:** Fusang + NJMerge-2 on 100 taxa, built from ten 10-taxon subsets. It was "slightly worse" than IQ-TREE and RAxML.
- **Zaharias, Grosshauser & Warnow**, *J Comput Biol* 2022, doi:10.1089/cmb.2021.0383.
  - DNN quartets + QMC were worse than MP, NJ and ML.
  - Their explanation: local-then-amalgamate loses the benefit of taxon sampling.
  - They predicted that DNN subsets plus supertree methods are "unlikely to be generally successful."
  - Our result supports that prediction at 50-taxon subsets with GTM.
- **Not found:** no work combining Phyloformer, Phyloformer 2 or NeuralNJ with GTM, NJMerge, TreeMerge or a supertree method, and no de novo DL tree estimator evaluated at 1000 taxa.
- **Course readings:**
  - Sapoval et al., *Nat Commun* 2022, doi:10.1038/s41467-022-29268-7.
  - Braichenko, Borges & Kosiol, *GBE* 2025, doi:10.1093/gbe/evaf177.
  - Tang et al., *Bioinformatics Advances* 2024, doi:10.1093/bioadv/vbae022.
  - Zhang et al., *MBE* 2025 = NeuralNJ.

## 3. Setup

**Data:** MAGUS-paper datasets (Illinois Data Bank IDB-2643961) from `cs581/code/setup.sh`, true alignments, binary true trees (`rose.tt`, `true_tree.tre`).
- Conditions: ROSE 1000M2, 1000L1, 1000S3, and RNASim 1000.
- Replicates: **R0–R2 for ROSE and R0–R1 for RNASim1000** (5 were planned; we stopped at the time budget).
- Phyloformer in the pipeline: **R0 only** for the three ROSE conditions, because of its CPU cost.

**Pipeline** (`code/run_rep.py`, one replicate, restartable, single-threaded; four replicates ran in parallel):
1. FastTree 2 (`-nt -gtr`) on the full alignment gives the guide tree, which is also the "FastTree full" baseline.
2. IQ-TREE 3.1.4, GTR+G, `--fast`, on the full alignment is the "IQ-TREE full" baseline.
3. Centroid-edge decomposition of the guide tree into subsets of **≤50 taxa** (same code as `cs581/gtm/code/sim.py`). That gave 27–30 subsets of 21–50 taxa (mean ≈ 34). We chose 50 because both DL models were trained on 50-taxon trees.
4. Subset trees from each method:
   - **FT**: FastTree GTR.
   - **IQ**: IQ-TREE GTR+G `--fast`.
   - **NJ**: FastME NJ on JC69 distances.
   - **BME**: FastME BME+NNI+SPR on JC69 distances; this is the classical counterpart of the Phyloformer tree search.
   - **NNJ**: NeuralNJ argmax.
   - **PF**: Phyloformer `pf.ckpt` distances → FastME BME+NNI+SPR. DNA letters A/C/G/T go to the amino-acid channels with the same letters; sites are split into windows of 256 and the distances averaged. `code/pfdist.py`.
5. Merge with GTM (github.com/vlasmirnov/GTM @18e3bc9, default mode), using the FastTree guide tree.
6. Baseline in the published setting (Park et al. 2021): GTM with IQ-TREE or FastTree subsets of **≤500** taxa, R0 only.

**Metric and statistics:**
- FN (missing-branch) rate against the true tree, using `phylo.py` from the GTM session (validated there against the published numbers).
- Comparisons are paired by replicate: mean difference, W/T/L, two-sided Wilcoxon signed-rank test.
- Every pipeline subset tree is also scored against the induced true tree, which is the small-tree validation on the true alignment.

**Other experiments:**
- **Estimated-alignment small-tree validation** (`code/val_small.py`): the 2 largest ≤50-taxon centroid subsets of the true tree, realigned from raw sequences with MAFFT L-INS-i. R0–R2, 4 conditions, 24 subsets.
- **Positive control** (`code/indomain_nnj.py`): NeuralNJ's own test set, 50 taxa at 256, 512 and 1024 sites, 10 alignments per length.
- **Domain-shift diagnostic** (`code/diag_nnj.py`): 5 subsets each from 1000M2 R0 and 1000S3 R0, in three versions:
  - "rose": the ROSE true alignment as is;
  - "nogap": every column containing a gap removed;
  - "alisim": IQ-TREE fits GTR+G4 and branch lengths on the fixed true subtree, then AliSim simulates 1000 sites with no indels.

## 4. Results

### 4.1 Full 1000-taxon trees (FN %, mean over replicates)

| method | 1000M2 (n=3) | 1000L1 (n=3) | 1000S3 (n=3) | RNASim1000 (n=2) |
|---|---|---|---|---|
| FastTree (full) | 9.37 | 9.33 | 9.82 | 14.24 |
| IQ-TREE `--fast` (full) | 9.53 | 10.17 | 9.82 | 15.30 |
| GTM + IQ subsets ≤500 (published setting; n=1) | 9.07 | 10.27 | 10.58 | 12.84 |
| GTM + FT subsets ≤500 (n=1) | 10.58 | 10.27 | 12.50 | 14.34 |
| GTM + IQ subsets ≤50 | 10.10 | 11.01 | 10.93 | 14.39 |
| GTM + FT subsets ≤50 | 10.91 | 12.62 | 11.88 | 14.64 |
| GTM + BME(JC) subsets ≤50 | 14.40 | 15.30 | 14.17 | 18.61 |
| GTM + NJ(JC) subsets ≤50 | 15.81 | 16.01 | 14.77 | 18.25 |
| **GTM + NeuralNJ subsets ≤50** | **18.06** | **19.66** | **18.71** | **23.52** |
| **GTM + Phyloformer subsets ≤50 (R0 only)** | **16.94** | **18.83** | **17.24** | – |

**Paired comparisons, pooled over conditions** (B − A in FN points; negative means B is better):

| A | B | n | mean diff | W/T/L (B vs A) | Wilcoxon p |
|---|---|---|---|---|---|
| GTM+IQ@50 | GTM+NeuralNJ@50 | 11 | +8.31 | 0/0/11 | 0.00098 |
| GTM+FT@50 | GTM+NeuralNJ@50 | 11 | +7.35 | 0/0/11 | 0.00098 |
| GTM+BME@50 | GTM+NeuralNJ@50 | 11 | +4.32 | 0/0/11 | 0.00098 |
| FastTree full | GTM+NeuralNJ@50 | 11 | +9.30 | 0/0/11 | 0.00098 |
| IQ-TREE full | GTM+NeuralNJ@50 | 11 | +8.84 | 0/0/11 | 0.00098 |
| GTM+IQ@50 | GTM+PF@50 | 3 | +6.78 | 0/0/3 | 0.25 |
| GTM+BME@50 | GTM+PF@50 | 3 | +3.02 | 0/0/3 | 0.25 |
| FastTree full | GTM+PF@50 | 3 | +7.25 | 0/0/3 | 0.25 |
| FastTree full | GTM+IQ@50 | 11 | +0.99 | 2/0/9 | 0.0049 |
| IQ-TREE full | GTM+IQ@50 | 11 | +0.52 | 3/0/8 | 0.083 |
| IQ-TREE full | GTM+IQ@500 | 4 | −0.98 | 3/0/1 | 0.25 |

The per-condition rows are in `results/SUMMARY.md`. Every condition goes the same way as the pooled rows.

**A side result that matters for any DL+DTM design: small subsets cost accuracy even with ML subset trees.**
- GTM with ≤50-taxon IQ-TREE subsets is about 1 point worse than FastTree on the full data (9/11 replicates worse, p = 0.005).
- With ≤500-taxon subsets (the published setting), GTM+IQ matches or beats full IQ-TREE (n = 1 per condition).
- So the subset size DL can handle today (≤50–100) already starts the merged tree at a disadvantage, before any DL error.

### 4.2 Small-tree validation: subset trees vs the induced true tree (FN %)

**True alignment** (all pipeline subsets of 21–50 taxa from R0–R2; n = 85 / 85 / 88 / 61 subsets, Phyloformer R0 only with n = 30 / 29 / 27):

| method | 1000M2 | 1000L1 | 1000S3 | RNASim1000 | s/subset (1 thread) |
|---|---|---|---|---|---|
| IQ (`--fast`) | 8.93 | 10.07 | 9.59 | 12.64 | 0.9 |
| FT | 9.83 | 11.68 | 10.50 | 13.05 | 1.5 |
| BME (JC) | 13.67 | 14.71 | 13.00 | 17.02 | <0.1 |
| NJ (JC) | 15.04 | 15.39 | 13.66 | 16.91 | <0.1 |
| NeuralNJ | 17.99 | 19.76 | 18.53 | 22.42 | 31 |
| Phyloformer + FastME (R0) | 15.96 | 18.68 | 15.88 | – | 50 |

Pooled over all subsets of R0–R2 (n = 319): NeuralNJ is worse than IQ by +9.3 points (W/T/L 26/25/268, p = 4e-44) and worse than NJ(JC) by +4.4 points (p = 1e-15). Phyloformer (R0, n = 86) is worse than IQ by +7.0 (p = 1e-9) and worse than NJ(JC) by +2.4 (p = 8e-4).

**MAFFT L-INS-i estimated alignment** (24 subsets, 47–50 taxa, R0–R2):

| method | 1000M2 | 1000L1 | 1000S3 | RNASim1000 | pooled |
|---|---|---|---|---|---|
| IQ | 8.00 | 11.22 | 8.56 | 9.90 | 9.42 |
| FT | 11.16 | 12.57 | 8.92 | 10.62 | 10.82 |
| BME (JC) | 15.64 | 17.18 | 11.48 | 15.01 | 14.82 |
| NJ (JC) | 16.53 | 18.26 | 10.74 | 14.64 | 15.04 |
| NeuralNJ | 20.38 | 21.97 | 15.22 | 22.32 | 19.97 |
| Phyloformer + FastME | 20.17 | 21.22 | 14.02 | 32.47 | 21.97 |

- NeuralNJ vs IQ: +10.6 points, 0/2/22, p = 4e-5.
- Phyloformer vs IQ: +12.6 points, 0/1/23, p = 3e-5.
- Phyloformer collapses on estimated RNASim alignments (up to 50% FN on one subset).

### 4.3 Positive control: NeuralNJ in its own training distribution (our CPU harness)

NeuralNJ's released test set (GTR+I+G AliSim, 50 taxa), 10 alignments per length:

| sites | NeuralNJ | IQ-TREE `--fast` | FastTree | BME (JC) | NJ (JC) | NeuralNJ vs IQ (W/T/L, p) |
|---|---|---|---|---|---|---|
| 256 | **29.8** | 37.9 | 45.5 | 54.0 | 57.4 | 6/2/2, p = 0.38 |
| 512 | **15.5** | 17.4 | 20.6 | 34.9 | 37.2 | 5/0/5, p = 0.76 |
| 1024 | 8.5 | **7.0** | 10.9 | 26.2 | 27.7 | 4/1/5, p = 0.47 |

- In distribution, NeuralNJ is competitive with ML and much better than distance methods, which is consistent with the paper.
- So the wrapper is not broken. The 2–4× degradation relative to IQ-TREE on ROSE and RNASim is out-of-distribution failure.

### 4.4 Diagnostic: is it indels, alignment length, or the trees?

Mean FN (%) over 5 subsets each of 1000M2 R0 and 1000S3 R0 (≈1000 gap-stripped columns per subset):

| data version | 1000M2: IQ / FT / BME / NeuralNJ | 1000S3: IQ / FT / BME / NeuralNJ |
|---|---|---|
| ROSE true alignment | 7.6 / 7.6 / 12.8 / **17.7** | 10.4 / 12.2 / 14.9 / **14.9** |
| gap-containing columns removed | 6.5 / 6.9 / 11.7 / **17.8** | 10.4 / 11.8 / 16.7 / **12.6** |
| AliSim GTR+G re-simulation on the true subtree, 1000 sites, no indels | 9.1 / 11.5 / 13.8 / **19.0** | 13.0 / 15.2 / 20.2 / **18.0** |

- Removing indel signal or re-simulating under a GTR model *without indels*, at an alignment length inside NeuralNJ's training range, does not rescue NeuralNJ on the 1000M2 subtrees.
- The remaining shift is in the **trees**: topology shape and branch lengths of centroid subsets of a 1000-taxon ROSE tree, as opposed to NeuralNJ's simulated 50-taxon training trees.
- On 1000S3 the gap is smaller.
- Mechanism (not tested): a centroid subset is a clade of a large tree, with a long stem, deep internal structure and a particular distribution of short internal branches. Its branch-length distribution likely differs from the training trees.

### 4.5 Runtime (seconds per replicate, 1 thread, with 4–5 jobs sharing 4 cores)

| method | 1000M2 | 1000L1 | 1000S3 | RNASim1000 |
|---|---|---|---|---|
| FastTree full | 128 | 145 | 87 | 238 |
| IQ-TREE `--fast` full | 168 | 235 | 137 | 177 |
| GTM+IQ@50 (guide + subsets + merge) | 150 | 167 | 129 | 259 |
| GTM+NeuralNJ@50 | 842 | 897 | 679 | 2029 |
| GTM+Phyloformer@50 | 1496 | 1566 | 1462 | – |

- On CPU, the DL subset step is 5–10× slower than the entire FastTree run, and 30–50 s per 50-taxon subset.
- NeuralNJ's cost grows with alignment length: about 70 s per subset on RNASim, where gap-stripped subsets have ~2000+ columns.
- Phyloformer on one thread is nearly as fast as on 4 threads (6.9 vs 8.6 s per 40×256 window). bf16 autocast gives NaNs.
- The GTM merge itself takes 0.2–0.5 s.

## 5. Caveats

- **Replicate counts.** Fewer replicates than planned (3/3/3/2, Phyloformer 1/1/1/0). The DL deficit is 4–10 points with 0 wins in 11 paired replicates and hundreds of paired subsets, so more replicates would not change the sign.
- **IQ-TREE `--fast`** was used everywhere for CPU budget. Full IQ-TREE would widen the gap to DL, not close it.
- **Phyloformer** was used out of domain by construction: a protein model fed DNA. Its result says nothing about a nucleotide-trained Phyloformer.
- **NeuralNJ** ran only in argmax mode, without RL fine-tuning or MC search; the paper's NeuralNJ-RL and NeuralNJ-MC variants need the raxml binding. Per the paper those variants are more accurate, but much slower and still on the same learned prior.
- **Subset-level p-values** treat subsets from one replicate as independent, which they are not quite. The replicate-level tests in §4.1 do not have this problem.
- **Wall-clock times** come from a shared, loaded machine. Compare within a column, not as absolute numbers.

## 6. What a 4-week project could be, and the risks

**Option A, recommended if the student wants DL: "Why DL tree estimators fail on DTM subsets, and does fine-tuning fix it?"**
- **Week 1:** the pipeline and harness exist (this directory). Extend to ≥10 replicates and the published GTM data (RNASim1000, Cox1-HET).
- **Week 2:** characterize the shift. Compare branch-length and tree-shape statistics (stemminess, tree height, internal/pendant ratio) of centroid subsets against NeuralNJ's training trees.
- **Week 3:** fine-tune NeuralNJ, or a small nucleotide Phyloformer, on CPU, using AliSim data simulated on *centroid subsets of simulated 1000-taxon trees*. Training labels are free.
- **Week 4:** rerun DL+GTM.
- **Success criterion:** GTM+DL@50 ≤ GTM+IQ@50.
- **Chance:** about 15–25% to reach parity with IQ-TREE subsets, under 10% to beat full-data FastTree or IQ-TREE.
- **Main risks:**
  1. CPU training: Phyloformer was trained on 6× A100s, and NeuralNJ's training config assumes a GPU. A CPU fine-tune of a 300k–1.3M-parameter model on maybe 10⁴ examples is feasible but may be too little.
  2. Even a perfect fine-tune is bounded by the small-subset penalty in §4.1: GTM+IQ@50 is already +1 point worse than FastTree full.
  3. NeuralNJ's raxml binding does not build on current toolchains (we used a shim).

**Option B, a safe negative-result paper:** the evaluation in §4, extended to Fusang (≤40 taxa, CPU, MIT) and to NJMerge/TreeMerge alongside GTM. Framed as a test of Zaharias et al. 2022's prediction at 50-taxon subsets. The result is near-certain, but the novelty is modest.

**Not recommended:** expecting pretrained DL + GTM to beat FastTree or IQ-TREE at 1000 taxa. The evidence here says no, by 7–10 FN points.

## 7. Reproduce

```bash
bash cs581/code/setup.sh            # FastTree, IQ-TREE 3, data
# Phyloformer: git clone https://github.com/lucanest/Phyloformer /opt/src/Phyloformer (commit 2cdd782); pip install torch (CPU)
# NeuralNJ:    git clone https://github.com/ZhangXinru99/NeuralNJ /opt/src/NeuralNJ (commit 3b59adb)
#              micromamba env /opt/mm/root/envs/nnj: python=3.10 + torch(cpu) "numpy<2" ete3 dendropy fvcore einops biopython tensorboard
#              raxmlpy shim: /opt/src/raxmlpy_shim/raxmlpy = pure-Python part of RAxMLpy/raxmlpy/core.py (treestr_to_tuples, utree2rtree_guided)
# GTM:         git clone https://github.com/vlasmirnov/GTM /opt/gtm/GTM_src && git -C /opt/gtm/GTM_src checkout 18e3bc9
cd cs581/dldtm/code
python3 run_rep.py 1000M2 R0 /opt/dldtm_runs 50 IQfull,FT,IQ,BME,NJ,NNJ,PF   # one replicate
bash run_queue.sh /opt/dldtm_runs 50 IQfull,FT,IQ,BME,NJ,NNJ "R1 R2 R3 R4"   # 4 at a time
bash run_queue.sh /opt/dldtm_runs 500 FT,IQ "R0"                            # published-setting baseline (run after the 50 queue for the same reps)
python3 val_small.py 1000M2 R0 /opt/dldtm_runs/val_est FT,IQ,NJ,BME,NNJ,PF   # estimated-alignment validation
python3 indomain_nnj.py /opt/dldtm_runs/indomain 10 ../results/indomain_nnj.jsonl
python3 diag_nnj.py 1000M2 R0 /opt/dldtm_runs/diag 5
python3 summarize.py /opt/dldtm_runs ../results
```

Note: do not run the ≤50 and ≤500 queues on the same replicate at the same time, because both write the same full-tree files. `run_queue.sh` throttles on the number of running `run_rep.py` processes but does not lock replicates.
