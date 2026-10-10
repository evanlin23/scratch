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
