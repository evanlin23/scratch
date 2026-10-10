# Rogue taxa and MSA / tree accuracy: a pilot (CS581 project idea)

## 1. Question and source

From the instructor's project list (`cs581/notes/CS581-project-suggestions.txt`, lines 59-61):

> Evaluate methods for tree construction and alignment estimation in the presence of rogue taxa
> (that are nevertheless homologs); how does the inclusion of the rogue taxa impact the accuracy
> of the alignment or tree restricted to the non-rogue taxa? Also consider whether removing rogue
> taxa before alignment and tree estimation improves accuracy.

This pilot asks whether the idea can carry a 4-week project. It measures (A) the damage that
rogues do to the non-rogue taxa, and (B) whether a cheap reference-free filter, followed by
UPP-style HMM add-back, removes that damage.

## 2. Prior art (short version; full list with DOIs in `results/prior_art.md`)

**Rogue identification is tree-centric.**
- RogueNaRok (Aberer, Krompass & Stamatakis 2013, *Syst Biol*, doi:10.1093/sysbio/sys078) and its RBIC predecessor (Pattengale et al. 2011, doi:10.1109/TCBB.2011.28).
- The leaf stability index (Thorley & Wilkinson 1999, doi:10.1006/jtbi.1999.0999).
- Safe taxonomic reduction and reduced consensus (Wilkinson 1994/1995/1996).
- Smith's information-theoretic `Rogue` (2022, *Syst Biol*, doi:10.1093/sysbio/syab099).
- TreeShrink (Mai & Mirarab 2018, *BMC Genomics*, doi:10.1186/s12864-018-4620-2), which looks for long branches that inflate the tree diameter.
- PhylteR (Comte et al. 2023, *MBE*, doi:10.1093/molbev/msad234) and its precursor Phylo-MCOA (de Vienne et al. 2012, doi:10.1093/molbev/msr317), which work on multi-gene distance matrices.

All of these are judged on tree or consensus quality. None is judged on the alignment.

**Long branches and taxon sampling.**
- LBA was first described by Felsenstein (1978, doi:10.2307/2412923); Bergsten (2005, doi:10.1111/j.1096-0031.2005.00059.x) reviews it.
- Adding taxa that break up long branches usually *helps*: Hillis 1996 (doi:10.1038/383130a0), Zwickl & Hillis 2002 (doi:10.1080/10635150290102339), Pollock et al. 2002 (doi:10.1080/10635150290102357).
- This argues against reflexive removal of taxa.

**Problem sequences in MSA.**
- UPP (Nguyen et al. 2015, *Genome Biol*, doi:10.1186/s13059-015-0688-z) and UPP2 (Park et al. 2023, doi:10.1093/bioinformatics/btad007) keep fragmentary sequences out of the backbone and add them with HMMs.
- MAGUS+UPP (Smirnov & Warnow 2021, *Syst Biol*, doi:10.1093/sysbio/syaa058), WITCH (Shen et al. 2022, doi:10.1089/cmb.2021.0585) and EMMA (Shen et al. 2023, doi:10.1186/s13015-023-00247-x) take the same route.
- These methods filter by length (fragments), not by divergence or rogue-ness, and they score all sequences, not the retained subset.

**Closest measurements of the "core subset" effect.**
- Sievers et al. 2013 (doi:10.1093/bioinformatics/btt093) and Boyce et al. 2014 (doi:10.1073/pnas.1405628111) score a fixed reference subset as more homologs are added. The variable is how many homologs, not how divergent they are.
- OD-seq (Jehl et al. 2015, doi:10.1186/s12859-015-0702-1) finds that removing outliers barely changes benchmark accuracy.
- AmpliPhy (Kim, Gil, Katoh & Dessimoz, bioRxiv 2026, doi:10.64898/2026.01.26.701724) adds homologs and finds better gene trees but almost no change in the alignment of the original sequences.

**Cautions on filtering.**
- Tan et al. 2015 (doi:10.1093/sysbio/syv033) found that column filtering usually makes single-gene trees worse.
- TAPER (Zhang et al. 2021, doi:10.1111/2041-210X.13696) warns that divergence is not error.

**Gap (moderate confidence: web and Crossref search only).** No study I found measures SP error
*restricted to the non-rogue taxa* as a function of including long-branch homologs, or tests whether
removing them before alignment and adding them back afterwards recovers that accuracy.

## 3. Design

### Data

The base data are the MAGUS paper's ROSE **1000M4** replicates: 1000 DNA taxa, mean pairwise p-distance
about 0.50, true alignment about 2.6k columns. 1000M2 is near saturation (mean p about 0.69), so a "long branch"
there is barely distinguishable from the rest.

There are three rogue settings, all with k = 50 rogues (5%). Rogue lists are in `results/detect/*.json`.

| setting | definition | n | rogue divergence |
|---|---|---|---|
| `M4lb` | (a) natural long-branch homologs: the 50 ROSE leaves with the longest terminal branch in the true model tree; no new sequences | 10 reps (R0-R9) | min terminal length about 5x the median; median nearest-neighbour p 0.30 vs 0.10 for others |
| `M4inj0.5` | (b) injected: 50 extra leaves, each attached at a uniformly random node (leaf or internal) of the true ROSE model tree, evolved from that node's **true ROSE sequence** (`rose.aln.true.internal.fasta`) with AliSim (IQ-TREE 3.1.4; K80, kappa=2, +G4 alpha=1, indel rate 0.05/0.05, geometric indel length mean 3) along a pendant branch of 0.5 subs/site | 10 reps | nearest p about 0.33 |
| `M4inj1.5` | same, pendant branch 1.5 | 5 reps | nearest p about 0.52 |

In the injected settings, the ROSE leaves' true alignment is unchanged; this is checked by an
assertion in `code/simulate.py`. A rogue residue goes into its anchor's column when AliSim keeps it
homologous; a rogue insertion gets its own column. So the reference for the non-rogue taxa is exactly
the paper's reference.

### Methods

- **Aligners:** MAFFT 7.505 `--auto` (FFT-NS-2 at this size); MAGUS(Fast) with the paper's flags (`gcmx.e2e_bench.magus_flags(25)`); PASTA 1.8.3 with 3 iterations.
- **Trees:** FastTree 2 `-nt -gtr -nosupport`.
- **Scoring:** every score is on the non-rogue taxa C only. Alignment error is the average of SPFN and SPFP from FastSP, computed on the induced sub-alignments of the estimate and the truth. Tree error is the FN rate (with FP) against the true tree restricted to C.

### Experiment A conditions, per aligner

| condition | alignment | tree |
|---|---|---|
| `all` | aligner on all 1050 (or 1000) taxa | FastTree on it |
| `all-treeC` | same alignment, rogue rows deleted | FastTree on C (isolates the effect of the alignment) |
| `oracle` | aligner on C only (known rogues removed) | FastTree on C |
| `oracle+add` | oracle alignment, with rogues added back by HMM | FastTree on all |
| `true` / `true-C` | true alignment with or without rogues | FastTree (isolates the tree-estimation effect) |

### Experiment B: reference-free detectors

Each detector runs once per instance.

- `ts`: TreeShrink 1.4 (pure Python) at its default alpha = 0.05, on a FastTree tree of an initial MAFFT alignment.
- `pd`: mean p-distance to the 5 nearest neighbours in the initial MAFFT alignment; flag when the robust z (median/MAD) is above 3.5.
- `hmm`: UPP-style backbone of 100 random sequences within 25% of the median length, aligned with MAFFT-L-INS-i, then `hmmbuild`. Every sequence is scored with `hmmsearch --max`; flag when bit score per residue has robust z below -3.5.

**Pipeline:** detect → align all minus flagged → add the flagged sequences back with `hmmbuild --hand`
(each backbone column is a match state) and `hmmalign`, so every backbone homology is unchanged and
each insertion gets its own columns → FastTree on all taxa. The `D-drop` variant leaves the flagged
taxa out of the tree permanently. It is scored on C minus the flagged set, paired against `all`
restricted to the same taxa.

**Statistics:** results are paired by replicate. Each table gives the mean difference, W/T/L (A better/tied/worse) and a
two-sided Wilcoxon signed-rank p (scipy; exact for small n).

## 4. Results

All numbers are in percent. The full tables are in `results/SUMMARY.md` (made by `code/summarize.py` from
`results/rows_*.jsonl`). Compute ran out before the planned 10 replicates, so the replicate counts are:
- MAFFT: 7 (M4inj0.5), 7 (M4lb) and 3 (M4inj1.5).
- MAGUS: 2 replicates, M4inj0.5 only.
- PASTA: 1 replicate, M4lb only.

These are small pilot samples. With n = 7 the smallest possible two-sided Wilcoxon p is 0.016.

### 4.1 Experiment A: how much do rogues hurt the non-rogue taxa?

Alignment error is the average of SPFN and SPFP, on the non-rogue taxa only.

| setting | aligner | n | with rogues (`all`) | rogues removed (`oracle`) | diff | W/T/L | p | control: 50 random non-rogues removed |
|---|---|---|---|---|---|---|---|---|
| M4inj0.5 | MAFFT | 7 | 4.70 | 4.38 | **-0.31** | 5/0/2 | 0.16 | +0.03 (3/0/4) |
| M4inj1.5 | MAFFT | 3 | 4.22 | 4.16 | -0.07 | 2/0/1 | 0.75 | +0.04 |
| M4lb | MAFFT | 7 | 4.30 | 4.05 | -0.25 | 5/0/2 | 0.16 | **-0.25 (6/0/1)** |
| M4inj0.5 | MAGUS | 2 | 0.85 | 0.78 | -0.07 | 1/0/1 | - | not run |
| M4lb | PASTA | 1 | 0.84 | 0.90 | +0.06 | 0/0/1 | - | not run |

Tree error is the FastTree FN rate, again on the non-rogue taxa only.

| setting | aligner | n | `all` | `oracle` | diff (p) | rogue rows dropped before FastTree (`all-treeC`) | oracle + HMM add-back |
|---|---|---|---|---|---|---|---|
| M4inj0.5 | MAFFT | 7 | 7.75 | 6.92 | **-0.83 (6/0/1, p=0.031)** | 7.34 (-0.42, p=0.06) | 6.98 (-0.77, 6/1/0, p=0.031) |
| M4inj1.5 | MAFFT | 3 | 6.95 | 6.75 | -0.20 | 6.75 | 6.69 |
| M4lb | MAFFT | 7 | 6.32 | 6.26 | -0.06 (p=0.89) | 6.29 | 6.20 |
| M4inj0.5 | MAGUS | 2 | 5.92 | 6.02 | +0.10 | 5.87 | 5.92 |
| M4lb | PASTA | 1 | 4.96 | 4.54 | -0.42 | 5.07 | 4.75 |
| *any* | *true alignment* | 8 / 3 / 7 | 6.76 / 6.45 / 6.09 | 6.78 / 6.45 / 6.32 | ~0 (rogues even help slightly in M4lb, 5/2/0) | | |

What this shows:

1. **The damage to the alignment is small, and it depends on the aligner.**
   - MAFFT (FFT-NS-2) loses about 0.3 points (about 7% relative) on the non-rogues when 50 moderately divergent homologs are present (M4inj0.5). The effect is mostly in SPFN (-0.49 vs -0.13 for SPFP). It is not significant at n = 7.
   - Very divergent rogues (L = 1.5, nearest p about 0.5) do *less* harm than moderate ones. MAFFT seems to push them aside instead of letting them pull the others out of place.
   - MAGUS and PASTA are already below 1% error. With or without rogues they differ by less than 0.1 point in either direction, which is inside replicate noise. The MAGUS and PASTA samples are tiny (n = 2 and n = 1).
2. **Natural long-branch ROSE leaves (M4lb) have no effect beyond what removing *any* 50 taxa does.** Removing 50 random non-rogue taxa gives the same -0.25 improvement with MAFFT (6/0/1). So "fewer sequences" accounts for the M4lb oracle gain, not "fewer rogues".
3. **Trees get worse only through the alignment.**
   - With the *true* alignment, including the rogues does not change FastTree's FN on the non-rogues (6.76 vs 6.78). On M4lb it even helps slightly, which matches the taxon-sampling literature.
   - With MAFFT on M4inj0.5, the oracle gains 0.83 FN points (p = 0.031). About half of that gain comes from deleting the rogue rows before tree estimation (`all-treeC`). The other half comes from the cleaner alignment.
   - Realigning without the rogues and then adding them back by HMM keeps almost all of the gain (-0.77), even though the rogues stay in the tree.
   - No such effect shows up for MAGUS. PASTA's one replicate is mixed.

### 4.2 Experiment B: reference-free detection and filtering

Detection quality is averaged over replicates. k = 50 true rogues in every instance.

| setting | TreeShrink (alpha=0.05) flagged / prec / recall | p-dist robust-z flagged / prec / recall | HMM bitscore robust-z flagged / prec / recall |
|---|---|---|---|
| M4inj0.5 | 1.1 / 0.57 / 0.02 | 39.1 / 0.98 / 0.77 | **50.0 / 1.00 / 1.00** |
| M4inj1.5 | 3.7 / 0.67 / 0.07 | 50.7 / 0.99 / 1.00 | **50.0 / 1.00 / 1.00** |
| M4lb | 3.0 / 0.54 / 0.03 | 0.4 / 0.29 / 0.01 | 0.9 / 0.24 / 0.01 |

The next table shows the effect of each filter followed by HMM add-back, compared with no filtering, using MAFFT on M4inj0.5 (n = 7). Each detector's flagged set is removed before alignment and then added back by HMM.

| pipeline | aln error | vs `all` | tree FN | vs `all` (W/T/L, p) |
|---|---|---|---|---|
| no filter (`all`) | 4.70 | | 7.75 | |
| TreeShrink | 4.65 | -0.04 | 7.58 | -0.17 (3/4/0, 0.25) |
| p-distance | 4.46 | -0.24 | 7.09 | -0.66 (6/0/1, 0.031) |
| HMM score | **4.38** | **-0.31** | **6.98** | **-0.77 (6/1/0, 0.031)** |
| oracle (known rogues removed, then added back) | 4.38 | -0.31 | 6.98 | -0.77 |

On M4lb no detector flags more than about 3 taxa, so every filter has essentially no effect (|diff| ≤ 0.1). On M4inj1.5 every filter is within ±0.13 of no filter.

What this shows:
- TreeShrink at its default alpha flags almost nothing on a single 1000-taxon tree. It is built to compare branch lengths across many gene trees, and per-tree it is very conservative.
- The UPP-style HMM score finds injected rogues perfectly. For MAFFT on M4inj0.5 it recovers all of the oracle gain: the alignment matches the oracle exactly in all 7 replicates.
- The ROSE leaves with the longest branches do not stand out to any detector, and removing them does not help. Being on a long branch in the true tree does not by itself make a taxon a rogue for alignment.

### 4.3 Runtime (4 cores, wall-clock in seconds per run)

All lanes shared the machine, so MAGUS and PASTA ran under about 3x CPU oversubscription. Their times are upper bounds on idle-machine time; the published MAGUS log for 1000M2 shows 968 s.

| step | wall |
|---|---|
| MAFFT `--auto` (1000-1050 taxa) | 13-22 |
| FastTree `-nt -gtr` (1000 taxa) | about 120 |
| HMM add-back of 50 sequences (hmmbuild --hand + hmmalign) | 3-6 |
| p-distance detector (incl. initial MAFFT) | 19-30 |
| HMM detector (100-seq L-INS-i backbone + hmmsearch) | 55-80 |
| TreeShrink detector (incl. initial MAFFT + FastTree) | 100-210 |
| MAGUS(Fast), paper flags | 1790-2030 |
| PASTA 1.8.3, 3 iterations | 3390-3470 |

The full pilot used about 3.5 h of wall-clock: 17 complete MAFFT instances plus 1 partial; about 100 MAFFT alignments, including the initial alignments the detectors need; 4 MAGUS and 2 PASTA alignments; and about 205 FastTree runs.

## 5. Caveats

- **Small n:** MAGUS n = 2 and PASTA n = 1 show only that there is no large effect. They cannot rule out a small one.
- **The injected rogues may be easy to detect.**
  - They evolve under K80+G4 with uniform base frequencies on top of ROSE sequences, so their composition drifts away from that of the ROSE sequences. The perfect HMM detection may partly exploit this model mismatch.
  - A ROSE-native simulation would avoid the mismatch: rerun ROSE on the tree with extra long leaves (ROSE is not installed here).
- **The control is scored on a different set.** The random-removal control is scored on C minus its 50 removed taxa, not on all of C. It calibrates the effect of "align fewer taxa"; it is not an exact null for the oracle comparison.
- **One model condition, one rogue fraction:** only 1000M4 and 5% rogues. Higher rogue fractions, protein data and biological rogues (contaminants, paralogs, misannotated taxa) are untested.
- **A bug was found and fixed during the run.** Early TreeShrink output parsing read `output_summary.txt` and wrongly counted all candidates as flagged. Those rows were deleted and rerun with the fixed parser, and the flagged sets in `results/detect/` are the corrected ones.

## 6. Verdict: **unclear, leaning not promising as posed**

- **Alignment.** As posed ("does including rogues hurt the alignment of the non-rogues, and does removing them first help?"), the pilot mostly gives a null result for the methods that matter in this course.
  - MAGUS and PASTA change by less than 0.1 point.
  - MAFFT loses about 0.3 points to moderately divergent injected rogues.
  - Natural long-branch taxa do nothing that removing random taxa would not.
- **Trees.** The positive finding is narrower. With a fast progressive aligner (MAFFT FFT-NS-2), injected rogues raise FastTree FN on the non-rogues by about 0.8 points (p = 0.03, n = 7). A cheap UPP-style HMM filter with add-back removes that loss completely. With the true alignment the rogues do no harm to the tree. So the harm, where it exists, comes through the alignment, not the tree method.

### What a 4-week project could look like

To be worth doing, the project would have to find the regime where the effect is large:

1. **Week 1: a dose-response grid.** Use MAFFT, plus MAGUS on a subsample. Vary rogue fraction (5, 15, 30%), pendant length (0.25-2) and rogue type: long branch, fragment plus long branch, and a compositionally biased rogue. Re-simulate with ROSE or INDELible so the rogues share the base model.
2. **Week 2: realistic rogues on biological data.**
   - Inject real distant homologs, such as Rfam or CRW sequences from a sister family, into the 16S datasets.
   - Repeat with protein (BAliBASE / HomFam), where Sievers et al. saw core-set degradation.
3. **Week 3: methods.**
   - Compare UPP/WITCH/EMMA add-back against plain HMM add-back.
   - Compare iterated TreeShrink (realign, retree), RogueNaRok on bootstrap trees, and OD-seq.
   - Report detector ROC curves, not single thresholds.
4. **Week 4: write-up.** Use the paired-replicate design from this pilot, including the random-removal control. That control is what showed the M4lb "effect" to be spurious.

### Risks

- The most likely outcome for MAGUS, PASTA and WITCH is still "no effect": their divide-and-conquer or HMM structure already keeps rogues from distorting the other subsets.
- Compute is a constraint. A paired MAGUS or PASTA comparison at 1000 taxa costs about 1-2 h on this machine, so 10+ replicates per cell needs a cluster or smaller datasets.
- The answer depends on how rogues are simulated, so the simulation model would have to be defended.
- A positive result might hold only for weak aligners, which makes it less interesting.

## 7. Files

- `code/simulate.py`: builds the instances (inject / longbranch) and checks that the base true alignment is preserved.
- `code/tools.py`: aligner wrappers, FastTree, HMM add-back, and the three detectors.
- `code/run.py`: the restartable per-instance driver for Experiments A and B plus the random control.
- `code/trees.py`: Newick parsing, restriction, and FN/FP.
- `code/summarize.py`: builds `results/SUMMARY.md`.
- `code/lane.sh` and `code/collect.sh`: run the lanes and collect results.
- `results/rows_*.jsonl`: one row per (instance, aligner, condition).
- `results/detect/*.json`: flagged sets, rogue lists and instance info.
- `results/prior_art.md`: the literature notes with DOIs.

To reproduce: `bash cs581/code/setup.sh`, clone TreeShrink to `/opt/src/TreeShrink` (`pip install treeswift`),
then `python simulate.py ...` (see the docstring) and `./lane.sh mafft M4inj0.5_R0 ...`.
